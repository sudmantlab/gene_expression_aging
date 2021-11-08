import allel
import argparse
import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from glmnet import ElasticNet

from sklearn.exceptions import ConvergenceWarning
import warnings

from itertools import islice
from joblib import Parallel, delayed
import sys

from collections import defaultdict
from itertools import chain
from operator import methodcaller
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings("ignore", category=ConvergenceWarning)
def fit_clf(gene_gt, expr, clf, tr_indiv=None, t_indiv=None, indiv_age_map=None, sep=False):
    output = {}
    Xt = np.transpose(gene_gt[t_indiv].to_numpy())
    coef = []
    perf = [] # Contains r2 performance of age and genetic factors
    if indiv_age_map != None and not sep: # Add age as last predictor if age_reg
        X = np.concatenate((Xt, np.expand_dims(np.array([indiv_age_map[ind] for ind in tr_indiv]), axis=1)),1)
    if indiv_age_map != None and sep:
        ageX = np.expand_dims(np.array([indiv_age_map[ind] for ind in tr_indiv]), axis=1)
        ageclf = LinearRegression()
        ageclf.fit(ageX, np.transpose(expr[tr_indiv]))
        coef += [ageclf.coef_[0]]
        X = Xt

    clf.fit(X, np.transpose(expr[tr_indiv]))
    # calcuate degrees of freedom genetics 
    degrees = 0
    if indiv_age_map != None and not sep:
        # Remove age coef from df calculation
        Xa = Xt[:, np.array(clf.coef_)[:-1]!=0]
    else: # get predictors that have nonzero coefs (and calcuate df)
        Xa = Xt[:, np.array(clf.coef_)!=0]
    degrees = np.trace(np.matmul(np.matmul(Xa, np.linalg.inv(np.matmul(Xa.T,Xa)+np.identity(Xa.shape[1])*0.5*clf.lambda_best_)), Xa.T))

    coef = list(clf.coef_) + coef
    return None, clf.cv_mean_score_[clf.lambda_best_inx_], coef, degrees
    
def train_model(gt, row, covar, indiv, tr_old_indiv, tr_young_indiv=None, t_old_indiv=None, t_young_indiv=None, eval_other=False, indiv_age_map=None, reg_out=True, sep=False):
    output = {} # Running dictionary of all results (R^2, pvalues, model coefficients) for current gene
    gene_gt = gt.loc[(gt["Chrom"]==row["#chr"]) & 
                     (gt["Pos"]>=(row["start"]-5e5)) & 
                     (gt["Pos"]<=(row["start"]+5e5)), 
                     gt.columns.difference(["Pos", "Chrom"])] # indiv x positions
    if gene_gt.shape[0] == 0: # No variants to be tested, terminate early and results should be 0's
        print("Zero variants for gene %s" % row["gene_id"])
        output["invalid_gene"] = [row["gene_id"]]
        return output # return zero performance

    if reg_out: # Regress the covariates out from expression data
        skclf = LinearRegression()
        skclf.fit(np.transpose(covar[indiv].to_numpy()), np.transpose(row[indiv].to_numpy()))
        covar_r2 = skclf.score(np.transpose(covar[indiv].to_numpy()), np.transpose(row[indiv].to_numpy()))
        row[indiv] = row[indiv] - skclf.predict(np.transpose(covar[indiv].to_numpy()))
        output["covar_r2"] = [covar_r2]

    # Fit regularized linear model for gene expresison for all or old individuals
    clf = ElasticNet(alpha=0.5, n_splits=10, n_jobs=1)
    old_tr_r2, old_t_r2, old_coef, df = fit_clf(gene_gt, row, clf, tr_old_indiv, t_old_indiv, indiv_age_map=indiv_age_map, sep=sep)
    if indiv_age_map!=None: # evaluate r2 of gt and age
        pred = np.matmul(old_coef[:-1], gene_gt[t_old_indiv].to_numpy())
        pred = pred - np.mean(pred)
        gene_r2 = r2_score(row[t_old_indiv]-np.mean(row[t_old_indiv]), pred) # Variance explained or R^2 for genotype weights in lin model (after correcting for covariates) 
        age_r2 = r2_score(row[t_old_indiv]-np.mean(row[t_old_indiv]), np.array([indiv_age_map[ind] for ind in t_old_indiv])*old_coef[-1]) # Variance explained by age weights
        output["age_r2"] = [age_r2]; output["gene_r2"] = [gene_r2]

    if eval_other and not indiv_age_map!=None: # Evaluates r2 on other cohort
        output["young_perf_old_model"] = [clf.score(np.transpose(gene_gt[tr_young_indiv].to_numpy()), np.transpose(row[tr_young_indiv]))]

    # build model and evaluate on young cohort if exists
    if tr_young_indiv is not None:
        clf = ElasticNet(alpha=0.5, n_splits=10, n_jobs=1)
        young_tr_r2, young_t_r2, young_coef, df = fit_clf(gene_gt, row, clf, tr_young_indiv, t_young_indiv, indiv_age_map=indiv_age_map, sep=sep)
        if indiv_age_map is not None: # evaluate r2 of gt vs age
            pred = np.dot(gene_gt[t_young_indiv].to_numpy(), np.array(young_coef[:-1]).reshape(-1,1))
            pred = pred - np.mean(pred)
            gene_r2 = r2_score(row[t_young_indiv]-np.mean(row[t_young_indiv]), pred) # Variance explained or R^2 for genotype weights in lin model (after correcting for covariates) 
            age_r2 = r2_score(row[t_young_indiv]-np.mean(row[t_young_indiv]), np.array([indiv_age_map[ind] for ind in t_young_indiv])*young_coef[-1]) # Variance explained by age weights
            output["young_age_r2"] = [age_r2]; output["young_gene_r2"] = [gene_r2]
        if eval_other: # Evaluates r2 on other cohort
            output["old_perf_young_model"] = [clf.score(np.transpose(gene_gt[tr_old_indiv].to_numpy()), np.transpose(row[tr_old_indiv]))]
        output["young_t_r2"]=[young_t_r2]; output["young_coef"] = [young_coef]
        output["old_t_r2"]=[old_t_r2]; output["old_coef"] = [old_coef] 
    else: # only one cohort - full dataset and age as predictor
        output["t_r2"]=[old_t_r2]; output["coef"] = [old_coef]
    output["gene"] = [row["gene_id"]]
    return output

# Get percent variance explained by each part of model (covariates, genotype, age - if using)
def eval_model(gt, row, covar, indiv, tr_old_indiv, coef, tr_young_indiv=None, t_old_indiv=None, t_young_indiv=None, indiv_age_map=None, reg_out=True):
    output = {}
    gene_gt = gt.loc[(gt["Chrom"]==row["#chr"]) & 
                     (gt["Pos"]>=(row["start"]-5e5)) & 
                     (gt["Pos"]<=(row["start"]+5e5)), 
                     gt.columns.difference(["Pos", "Chrom"])] # indiv x positions
    if gene_gt.shape[0] == 0:
        print("Zero variants for gene %s" % row["gene_id"])
        output["covar_r2"]=output["gene_r2"]=output["age_r2"]=output["coef"]=[0]
        return output # return zero performance
    if reg_out: # Regress the covariates out from expression data
        skclf = LinearRegression()
        skclf.fit(np.transpose(covar[indiv].to_numpy()), np.transpose(row[indiv].to_numpy()))
        covar_r2 = skclf.score(np.transpose(covar[indiv].to_numpy()), np.transpose(row[indiv].to_numpy()))
        row[indiv] = row[indiv] - skclf.predict(np.transpose(covar[indiv].to_numpy()))

    gene_r2 = r2_score(row[tr_old_indiv], np.matmul(coef[:-1], gene_gt[tr_old_indiv].to_numpy())) # Variance explained or R^2 for genotype weights in lin model (after correcting for covariates) 
    if indiv_age_map != None:
        age_r2 = r2_score(row[tr_old_indiv], np.array([indiv_age_map[ind] for ind in tr_old_indiv])*coef[-1]) # Variance explained by age weights
        output["age_r2"] = [age_r2]
    output["covar_r2"]=[covar_r2]; output["gene_r2"] = [gene_r2]; output["coef"] = [coef]
    return output

if __name__=="__main__":
    gtex_dir = "/clusterfs/genomicdata/GTEx/eQTL_files/"
    data_dir = "/clusterfs/nilah/rkchung/data/expred/"

    parser = argparse.ArgumentParser()
    
    # add arguments
    parser.add_argument("--expr", "-i", type=str, help="Input expression file. If none provided, runs all tissues")
    parser.add_argument("--cohorts", "-c", action="store_true", help="Split individuals into old and young cohorts. Default: False, runs without age cohorts")
    parser.add_argument("--no_covar", "-nc", action="store_true", help="No covariate effects removed. Default: False, regress out covariates")
    parser.add_argument("--sep", "-s", action="store_true", help="Separate model for age and genetics. Default: False, joint model")
    parser.add_argument("--processed", "-p", action="store_true", help="Use regressed out covariates. Default: False, regress out covariates")
    parser.add_argument("--split_tt", "-tt", action="store_true", help="Splits individuals into train and test datasets")
    parser.add_argument("--run", "-r", type=str, help="Name of run for replicates.")
    # read arguments from the command line
    args = parser.parse_args()

    cohorts = args.cohorts # Split individuals into age cohorts or build one model with all individuals
    reg_out = not args.no_covar # Regress out covariates (PEER factors, genetic PCs, sex, batch info etc.)
    sep = args.sep # Separate models for age and genetic components
    processed = args.processed # Use processed PEER factors (age components removed) 
    split_tt = args.split_tt# training and testing split
    tissues = [args.expr]
    run = args.run # Run number (if repeating analysis) else None

    # Load genotype data for all individuals
    callset = allel.read_vcf(data_dir+"gtex.filt.vcf")
    gt = np.sum(callset["calldata/GT"], axis=-1).astype(np.int8)
    chrom = callset['variants/CHROM']
    pos = callset["variants/POS"]
    indiv = [name.split("_")[1] for name in callset['samples']]
    gt = pd.DataFrame(data=gt, columns=indiv)
    pos = pd.DataFrame(data=pos, columns=["Pos"])
    chrom = pd.DataFrame(data=chrom, columns=["Chrom"])
    gt = pd.concat((chrom, pos, gt), axis=1)
    gt["Chrom"] = gt["Chrom"].astype(str)

    if args.expr is None: # Provide name of tissue to analyze
        tissues = glob(gtex_dir+"GTEx_Analysis_v8_eQTL_expression_matrices/*.bed.gz")
    print(tissues)

    # Load age of each individual
    age_df = pd.read_csv(data_dir+"phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt", sep="\t")[["SUBJID", "AGE"]]
    age_df = age_df.rename(columns={"AGE": "age"})
    age_df = age_df.drop_duplicates()
    age_df["age_cohort"] = "young"
    age_df.loc[age_df["age"]>55, "age_cohort"] = "old"
    age_df["norm_age"] = (age_df["age"]-age_df["age"].mean())/age_df["age"].std()
    age_df.to_csv("ages.tmp")
    indiv_age_map = dict(zip(age_df.SUBJID, age_df.norm_age))

    # Split individuals into age cohorts
    old_indiv = list(age_df[age_df["age_cohort"]=="old"]["SUBJID"])
    young_indiv = list(age_df[age_df["age_cohort"]=="young"]["SUBJID"])

    old_r2 = []
    young_r2 = []
    start = 0

    print("Run num %s" % run)
    for i in range(len(tissues)):
        short = tissues[i].split("expression_matrices/")[1].split(".v8.normalized_expression")[0]
        output_name = data_dir+"output/%s%s%s%s%s%s%s" % ("res_", short, "_cohort" if cohorts else "_agereg", "_proc" if processed else "", "_sep" if sep else "", "_split" if split_tt else "", run if run is not None else "")
        print("Generating %s.npy" % output_name)
        if processed:
            adj_dir = "/global/scratch-old/ryo10244201/analysis/Predixcan_cov/"
            covar = pd.read_csv(adj_dir + short + "_covariates.csv", index_col=0)
            #expr = pd.read_csv(tissues[i], index_col=0)

            # remove age and transpose expr and covariate dataframes
            covar.index = pd.Index(covar["SUBJID"])
            covar.drop(["AGE", "SUBJID"], axis=1, inplace=True) 
            covar = covar.T
        else:
            covar = pd.read_csv(gtex_dir + "GTEx_Analysis_v8_eQTL_covariates/" + short + ".v8.covariates.txt", delim_whitespace=True)
            print("Starting Tissue %s" % short)
            
            # Get GTEx covariates and select top n PEER factors
            num_peer = 15
            remove_covar = ["InferredCov"+str(j) for j in range(num_peer+1, 61)]
            covar = covar[~covar.ID.isin(remove_covar)]

        # Read expression data
        expr = pd.read_csv(tissues[i], sep="\t")
        expr["#chr"] = expr["#chr"].str.split("chr").str[1]
        expr["gene_id"] = expr["gene_id"].str.split(".").str[0]

        indiv = list(set(expr.columns).intersection(set(indiv)))
        
        # All individual expression model
        if not cohorts:
            indiv = np.array(indiv)
            if split_tt:
                split = 0.2 
                perm = np.random.permutation(len(indiv))
                t_indiv = indiv[perm[:int(split*len(indiv))]]
                tr_indiv = indiv[perm[int(split*len(indiv)):]]
            else:
                perm = np.random.permutation(len(indiv))
                tr_indiv = indiv[perm]
                t_indiv = indiv[perm] # evaluate on training set

            # For each gene, train a PrediXcan-type gene expression model
            result_items = map(methodcaller('items'), Parallel(n_jobs=24)(delayed(train_model)(gt, row, covar, indiv, tr_indiv, t_old_indiv=t_indiv, indiv_age_map=indiv_age_map, reg_out=reg_out, sep=sep)
                    for g, row in tqdm(islice(expr.iterrows(), start, None), desc="Training Genes", total=(len(expr)-start))))
            results = defaultdict(list)
            for k, v in chain.from_iterable(result_items):
                results[k].extend(v)
            print("Average Heritability: %s+/-%s" % (np.mean(results["t_r2"]), np.std(results["t_r2"])))

        # Run cohort analysis
        else:  
            # Select young and old cohorts to match number of individuals in each
            old_indiv = np.array(list(set(indiv).intersection(set(old_indiv))))
            young_indiv = np.array(list(set(indiv).intersection(set(young_indiv))))
            cohort_perm = np.random.permutation(max(len(old_indiv), len(young_indiv)))
            if len(old_indiv) > len(young_indiv): # randomly subsample old cohort to match len of young cohort
                old_indiv = old_indiv[cohort_perm[:len(young_indiv)]]
            elif len(old_indiv) < len(young_indiv): # randomly subsample young cohort to match len of old cohort
                young_indiv = young_indiv[cohort_perm[:len(old_indiv)]]

            # Split individuals into train and test set
            if split_tt:
                split = 0.2 
                perm = np.random.permutation(len(old_indiv))
                t_old_indiv = old_indiv[perm[:int(split*len(old_indiv))]]
                tr_old_indiv = old_indiv[perm[int(split*len(old_indiv)):]]
                perm = np.random.permutation(len(young_indiv))
                t_young_indiv = young_indiv[perm[:int(split*len(young_indiv))]]
                tr_young_indiv = young_indiv[perm[int(split*len(young_indiv)):]]
            else:
                perm = np.random.permutation(len(old_indiv))
                tr_old_indiv = old_indiv[perm]
                perm = np.random.permutation(len(young_indiv))
                tr_young_indiv = young_indiv[perm]
                t_old_indiv = old_indiv[perm] # evaluate on training set
                t_young_indiv = young_indiv[perm]
            
            # For each gene, train a PrediXcan gene expression model
            result_items = map(methodcaller('items'), Parallel(n_jobs=24)(delayed(train_model)(gt, row, covar, indiv, tr_old_indiv, tr_young_indiv, t_old_indiv=t_old_indiv, t_young_indiv=t_young_indiv, reg_out=reg_out)
                    for g, row in tqdm(islice(expr.iterrows(), start, None), desc="Training Genes", total=(len(expr)-start))))
            results = defaultdict(list)
            for k, v in chain.from_iterable(result_items):
                results[k].extend(v)
        print("Heritability of older indivs: %s+/-%s" % (np.mean(results["old_r2"]), np.std(results["old_r2"])))
        print("Heritability of younger indivs: %s+/-%s" % (np.mean(results["young_r2"]), np.std(results["young_r2"])))
        np.save(output_name, results)
