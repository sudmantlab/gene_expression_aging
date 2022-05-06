import numpy as np
import pandas as pd

import allel
from sklearn.decomposition import PCA

import os

def load_expr_covar(pivus, processed, short, gtex_dir):
    expr = pd.read_csv(short, sep="\t")
    expr["#chr"] = expr["#chr"].str.split("chr").str[1]

    short = short.split("expression_matrices/")[1].split(
            ".v8.normalized_expression")[0]
    # combine PEER factors, PCA, sex into covar df
    if processed:
        adj_dir = os.getcwd()+"/../PEER/"
        #adj_dir = "/global/scratch/users/ryo10244201/analysis/Predixcan_cov/"
        covar = pd.read_csv(adj_dir + short + "_covariates.csv", index_col=0)

        # remove age and transpose expr and covariate dataframes
        covar.index = pd.Index(covar["SUBJID"])
        covar.drop(["AGE", "SUBJID"], axis=1, inplace=True) 
        covar = covar.rename(
                {"V6":"pcr", "V7":"platform", "V8":"sex"}, axis=1)
        covar = covar.T
    else:
        covar = pd.read_csv(gtex_dir + "GTEx_Analysis_v8_eQTL_covariates/" 
                            + short + ".v8.covariates.txt",
                            delim_whitespace=True)
        print("Starting Tissue %s" % short)
        
        # Get GTEx covariates and select top n PEER factors
        num_peer = 15
        remove_covar = ["InferredCov"+str(j) for j in range(num_peer+1, 61)]
        
        covar = covar[~covar.ID.isin(remove_covar)]
    return expr, covar, short

def load_data(data_dir, gtex_dir, age_reg_out=False):
    callset = allel.read_vcf(data_dir+"gtex.filt.vcf")
    if not age_reg_out:
        gt = np.sum(callset["calldata/GT"], axis=-1).astype(np.int8)
        chrom = callset['variants/CHROM']
        pos = callset["variants/POS"]
        indiv = [name.split("_")[1] for name in callset['samples']]
        gt = pd.DataFrame(data=gt, columns=indiv)
        pos = pd.DataFrame(data=pos, columns=["Pos"])
        chrom = pd.DataFrame(data=chrom, columns=["Chrom"])
        gt = pd.concat((chrom, pos, gt), axis=1)
        gt["Chrom"] = gt["Chrom"].astype(str)

    # Age of each individual
    age_df = pd.read_csv(data_dir+"phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt", sep="\t")[["SUBJID", "AGE"]]
    age_df = age_df.rename(columns={"AGE": "age"})
    age_df = age_df.drop_duplicates()
    age_df["age_cohort"] = "young"
    age_df.loc[age_df["age"]>55, "age_cohort"] = "old"
    age_df["norm_age"] = (age_df["age"]-age_df["age"].mean())/age_df["age"].std()
    indiv_age_map = dict(zip(age_df.SUBJID, age_df.norm_age))

    # Individuals in each age cohort
    #age_df = pd.read_csv("/clusterfs/nilah/rkchung/data/expred/Sample_age.csv", sep=";")
    old_indiv = list(age_df[age_df["age_cohort"]=="old"]["SUBJID"])
    young_indiv = list(age_df[age_df["age_cohort"]=="young"]["SUBJID"])
    return gt, indiv, old_indiv, young_indiv, indiv_age_map 

def load_pivus_data(pivus_dir, data_dir):
    # Load mapping of individuals in age 70 cohort to age 80 cohort
    map_df = pd.read_csv(pivus_dir+"Map_GlobalID.txt", sep="\t")
    map_df = map_df[["RNAseqID70", "RNAseqID80"]].sort_values("RNAseqID70").dropna()
    id_map = dict(zip(map_df.RNAseqID70, map_df.RNAseqID80))

    # Load sex of individuals
    sex_df = pd.read_csv(pivus_dir+"expCovariatesFinal.txt", sep="\t")[["Sex"]]

    # Load precomputed PEER factors
    peer = pd.read_csv(pivus_dir+"Peer_PIVUS.tsv", sep="\t")[["V"+str(i) for i in range(27,42)]]

    # Load individual genotype 
    gt = pd.read_csv(pivus_dir+"pilot70_impute2_af.05_info.4.hwe.meangens", sep="\t")
    gt[["Chrom", "Pos", "_", "_"]] = gt.SNP.str.split("_", expand=True)
    gt["Pos"] = gt["Pos"].astype(int)

    # Matched individuals in study
    young_indiv = sorted(list(set(map_df.RNAseqID70).intersection(gt.columns)))
    old_indiv = [id_map[iv] for iv in young_indiv]
    indiv = sorted(list(set(young_indiv).union(set(old_indiv))))

    ncomp = 5
    pca = PCA(n_components=ncomp)
    pccomp = pca.fit_transform(gt[sorted(list(set(gt.columns).intersection(indiv)))].T.to_numpy())
    pca_df = pd.DataFrame(pccomp, columns=["PC"+str(i) for i in range(ncomp)])
    pca_df = pd.concat((pca_df, pca_df), axis=0)
    pca_df.index = young_indiv+old_indiv
    pca_df = pca_df.sort_index()

    gt_old = gt[young_indiv].copy()
    gt_old.columns = old_indiv
    gt = pd.concat((gt, gt_old), axis=1)

    peer.index = pca_df.index
    covar = pd.concat((pca_df, peer, sex_df), axis=1).dropna()

    # Read in PIVUS normalized expression
    expr = pd.read_csv(pivus_dir+"processed_normcounts.txt", delim_whitespace=True)
    expr = expr.reset_index().rename(columns={'index': 'gene_id'})
    expr = expr.loc[expr["gene_id"].apply(lambda g: len(g.split("ENSG"))==2)]

    gene_tss = pd.read_csv(data_dir+"hg19_gene_tss.tsv", 
                           sep="\t"
                          ).rename({"Chromosome/scaffold name":"#chr", 
                                    "Gene stable ID version":"gene_id", 
                                    "Transcription start site (TSS)":"start"}, 
                          axis=1)
    gene_tss = gene_tss.drop_duplicates("gene_id")
    #gene_tss = gene_tss[["#chr", "start", "gene_id"]].drop_duplicates("gene_id")
    expr = expr.merge(gene_tss, on="gene_id")
    expr = expr[pd.to_numeric(expr['#chr'], errors='coerce').notnull()]


    # Filter by genes with GTEx read filter (>6 reads in at least 20% of samples)
    # and and TPM filter (>0.1 TPM in at least 20% of samples)
    gene_df = pd.read_csv(data_dir+"pivus_filtered_genes.txt")
    expr = expr.merge(gene_df, on="gene_id")

    # Make sure individuals match expression columns
    young_indiv = list(set(young_indiv).intersection(expr.columns))
    old_indiv = list(set(old_indiv).intersection(expr.columns))
    indiv = list(set(indiv).intersection(expr.columns))
    return gt, expr, covar.T, indiv, old_indiv, young_indiv
