import pandas as pd
import numpy as np
from glob import glob
import argparse

from sklearn.exceptions import ConvergenceWarning
import warnings

from itertools import islice
from itertools import chain

from tqdm import tqdm
from joblib import Parallel, delayed

from collections import defaultdict
from operator import methodcaller
from scipy import stats
import os

from utils.helpers import Seed
from utils.dataloader import load_data, load_pivus_data, load_covar
from utils.train import *


warnings.filterwarnings("ignore", category=ConvergenceWarning)

if __name__=="__main__":
    gtex_dir = "/clusterfs/genomicdata/GTEx/eQTL_files/"
    data_dir = "/clusterfs/nilah/rkchung/data/expred/"

    parser = argparse.ArgumentParser()
    metaseed = 423671
    
    # add arguments
    parser.add_argument("--expr", "-i", type=str, 
            help="Input expression file. If none provided, runs all tissues")
    parser.add_argument("--cohorts", "-c", action="store_true", 
            help="Split individuals into old and young cohorts. \
                    Default: False, runs without age cohorts")
    parser.add_argument("--no_covar", "-nc", action="store_true", 
            help="No covariate effects removed. \
                    Default: False, regress out covariates")
    parser.add_argument("--sep", "-s", action="store_true", 
            help="Separate model for age and genetics. \
                    Default: False, joint model")
    parser.add_argument("--sex", "-sex", action="store_true", 
            help="Include interaction term for age and sex \
                    Default: False")
    parser.add_argument("--pivus", "-pv", action="store_true", 
            help="Run on PIVUS dataset instead of GTEx. Default: False")
    parser.add_argument("--processed", "-p", action="store_true", 
            help="Use regressed out covariates. \
                    Default: False, regress out covariates")
    parser.add_argument("--split_tt", "-tt", action="store_true", 
            help="Splits individuals into train and test datasets")
    parser.add_argument("--run", "-r", type=str, 
            help="Name of run for replicates.")
    parser.add_argument("--bootstrap", "-b", type=str, 
            help="Number of bootstrap samples.")
    parser.add_argument("--saveinterm", "-si", action="store_true", 
            help="Save intermediate results when genes finished. \
                    Default: False, does not save intermediate results")
    # read arguments from the command line
    args = parser.parse_args()

    cohorts = args.cohorts # Split individuals into age cohorts or build one model with all individuals
    reg_out = not args.no_covar # Regress out covariates (PEER factors, genetic PCs, sex, batch info etc.)
    sep = args.sep # Separate models for age and genetic components
    sex = args.sex # Age interaction term for sex
    pivus = args.pivus # Age interaction term for sex
    bootstrap = int(args.bootstrap) if args.bootstrap != None else 1  # Number of bootstrap samples (none if no bootstrapping) 
                    # and provide stats about number of times beta nonzero coef
    processed = args.processed # Use processed PEER factors (age components removed) 
    saveinterm = args.saveinterm
    split_tt = args.split_tt# training and testing split
    tissues = [args.expr]
    age_reg_out = False
    run = args.run # Run number (if repeating analysis) else None

    if pivus:
        pivus_dir = "/global/scratch/users/ryo10244201/PIVUS/"
        gt, indiv, old_indiv, young_indiv = \
                                      load_pivus_data(pivus_dir, data_dir)
        tissues = ["pivus"]
    else:
        gt, indiv, old_indiv, young_indiv, indiv_age_map = \
                                                    load_data(data_dir, gtex_dir,
                                                              age_reg_out)
        if args.expr is None: # Provide name of tissue to analyze
            tissues = glob(gtex_dir+"GTEx_Analysis_v8_eQTL_expression_matrices/*.bed.gz")
        print(tissues)

    old_r2 = []
    young_r2 = []
    start = 0

    print("Run num %s" % run)
    for i in range(len(tissues)):
        short = tissues[i]

        if not pivus:
            expr, covar = load_expr_covar(pivus, processed, short, gtex_dir)

        if not age_reg_out:
            indiv = list(set(expr.columns).intersection(set(indiv)))
        else:
            indiv = list(expr.columns)[4:]
            print(indiv)
        expr["gene_id"] = expr["gene_id"].str.split(".").str[0]

        outputname = data_dir+"output/%s%s%s%s%s%s%s%s%s" % \
                ("res_", short,
                        "_agereg" if not cohorts else "_cohort", 
                        "_regoutage" if age_reg_out else "", 
                        "_proc" if processed else "", 
                        "_sep" if sep else "", 
                        "_split" if split_tt else "", 
                        "_sex" if sex else "", 
                        run if run is not None else "")
        print(outputname)

        seed = Seed(short, metaseed)

        if age_reg_out:
            model_res = np.load(
                    data_dir+"output/res_%s%s.npy" % 
                            (short, "_agereg" if not cohorts else ""), 
                    allow_pickle=True).item()

            # For each gene, evaluate PrediXcan gene expression model
            start = len(expr)-2
            result_items = map(
                    methodcaller('items'), 
                    Parallel(n_jobs=20)(
                        delayed(reg_out_age)(
                            row, indiv, indiv, indiv_age_map=indiv_age_map)
                        for g, row in tqdm(
                            islice(expr.iterrows(), start, None), 
                            desc="Eval Genes", total=(len(expr)-start))))
            results = defaultdict(list)
            for k, v in chain.from_iterable(result_items):
                results[k].extend(v)
            procexpr = expr.copy()
            procexpr.loc[start:, indiv] = np.array(results["procexpr"])
            del results["procexpr"]
            print("Average Var Explained: Age %s" % np.std(results["regoutage_r2"]))
            procexpr.to_csv(data_dir+"output/expr/"+short+".tsv", 
                            sep="\t", index=False)
        else:
            if not cohorts:
                indiv = np.array(indiv)
                if split_tt:
                    split = 0.2 
                    np.random.seed(seed.get(t="permute"))
                    perm = np.random.permutation(len(indiv))
                    t_indiv = indiv[perm[:int(split*len(indiv))]]
                    tr_indiv = indiv[perm[int(split*len(indiv)):]]
                else:
                    np.random.seed(seed.get(t="permute"))
                    perm = np.random.permutation(len(indiv))
                    tr_indiv = indiv[perm]
                    t_indiv = indiv[perm] # evaluate on training set
                # For each gene, train a PrediXcan gene expression model
                intermresdir = outputname + "_interm" if saveinterm else None
                if intermresdir!=None and not os.path.isdir(intermresdir):
                    os.mkdir(intermresdir)

                result_items = map(
                        methodcaller('items'), 
                        Parallel(n_jobs=24)(
                            delayed(train_model)( 
                                gt, row, covar, indiv, tr_indiv, 
                                t_old_indiv=t_indiv,indiv_age_map=indiv_age_map, 
                                reg_out=reg_out, sep=sep, bootstrap=bootstrap,
                                intermresdir=intermresdir, sex=sex, 
                                seed=seed)
                            for g, row in tqdm(
                                islice(expr.iterrows(), start, None),
                                desc="Training Genes", total=(len(expr)-start)))
                        )
                results = defaultdict(list)
                for k, v in chain.from_iterable(result_items):
                    results[k].extend(v)
                print("Average Heritability: %s+/-%s" % 
                        (np.mean(results["t_r2"]), np.std(results["t_r2"])))
            else:
                # Select young and old cohorts to match number of individuals in each
                old_indiv = np.array(list(set(indiv).intersection(set(old_indiv))))
                young_indiv = np.array(list(set(indiv).intersection(set(young_indiv))))

                cohort_perm = np.random.permutation(max(len(old_indiv), len(young_indiv)))

                # randomly subsample old cohort to match len of young cohort
                if len(old_indiv) > len(young_indiv): 
                    old_indiv = old_indiv[cohort_perm[:len(young_indiv)]]
                # randomly subsample young cohort to match len of old cohort
                elif len(old_indiv) < len(young_indiv): 
                    young_indiv = young_indiv[cohort_perm[:len(old_indiv)]]
                print(len(old_indiv))
                print(len(young_indiv))

                # Split individuals into train and test set
                if split_tt:
                    split = 0.2 
                    np.random.seed(seed.get(t="permute", old=True))
                    perm = np.random.permutation(len(old_indiv))
                    t_old_indiv = old_indiv[perm[:int(split*len(old_indiv))]]
                    tr_old_indiv = old_indiv[perm[int(split*len(old_indiv)):]]
                    np.random.seed(seed.get(t="permute", old=False))
                    perm = np.random.permutation(len(young_indiv))
                    t_young_indiv = young_indiv[perm[:int(split*len(young_indiv))]]
                    tr_young_indiv = young_indiv[perm[int(split*len(young_indiv)):]]
                else:
                    np.random.seed(seed.get(t="permute", old=True))
                    perm = np.random.permutation(len(old_indiv))
                    tr_old_indiv = old_indiv[perm]
                    np.random.seed(seed.get(t="permute", old=False))
                    perm = np.random.permutation(len(young_indiv))
                    tr_young_indiv = young_indiv[perm]
                    t_old_indiv = old_indiv[perm] # evaluate on training set
                    t_young_indiv = young_indiv[perm]
                print(expr.shape)
                
                # For each gene, train a PrediXcan gene expression model
                intermresdir = outputname + "_interm" if saveinterm else None
                if intermresdir!=None and not os.path.isdir(intermresdir):
                    os.mkdir(intermresdir)
                result_items = map(
                        methodcaller('items'), 
                        Parallel(n_jobs=24)(
                            delayed(train_model)(
                                gt, row, covar, indiv, tr_old_indiv, 
                                tr_young_indiv, t_old_indiv=t_old_indiv, 
                                t_young_indiv=t_young_indiv, reg_out=reg_out, 
                                bootstrap=bootstrap, intermresdir=intermresdir,
                                sex=sex, seed=seed)
                            for g, row in tqdm(
                                islice(expr.iterrows(), start, None), 
                                desc="Training Genes", total=(len(expr)-start))))
                # Compile all results from saved intermediate directory
                if intermresdir!= None:
                    result_items = []
                    files = glob(intermresdir+"/*")
                    for f in files:
                        res = dict(np.load(f, allow_pickle=True).item())
                        result_items += [res.items()] 

                results = defaultdict(list)
                for k, v in chain.from_iterable(result_items):
                    results[k].extend(v)

            print("Heritability of older indivs: %s+/-%s" % 
                    (np.mean(results["old_r2"]), np.std(results["old_r2"])))
            print("Heritability of younger indivs: %s+/-%s" %  
                    (np.mean(results["young_r2"]), np.std(results["young_r2"])))
        np.save(outputname, results)
