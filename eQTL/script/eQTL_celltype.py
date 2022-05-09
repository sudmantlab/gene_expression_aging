import os
from numpy.core.defchararray import join
import pandas as pd
from scipy.stats import entropy, ttest_ind, pearsonr, spearmanr, ttest_1samp
import numpy as np
from scipy.spatial import distance as dist
import time
import argparse
from JSD_utils import *
import statsmodels.formula.api as smf
from pysam import VariantFile, TabixFile
import biomart

parser = argparse.ArgumentParser(description='Snakemake script for calculating adjacent gene Spearman correlation.')
parser.add_argument('-t', '--test', action='store_true', 
                    default=False,
                    help='testing script with single Blood input')

parser.add_argument('-f', '--frame', action='store_true', 
                    default=False,
                    help='saving frame with single Blood input')

args = parser.parse_args()

server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')         
mart = server.datasets['hsapiens_gene_ensembl'] 
attributes = ['hgnc_symbol', 'ensembl_gene_id']
response = mart.search({'attributes': attributes})
data = response.raw.data.decode('ascii')
dict = {}
dict_rev = {}                                                
for line in data.splitlines():                                              
    line = line.split('\t')                                                 
    # The entries are in the same order as in the `attributes` variable                                              
    gene_symbol = line[0]                                                   
    ensembl_gene = line[1]                                         
                                                                            
    # Some of these keys may be an empty string. If you want, you can 
    # avoid having a '' key in your dict by ensuring the attributes
    # have a nonzero length before adding them to the dict   
    dict[ensembl_gene] = gene_symbol
    dict_rev[gene_symbol] = ensembl_gene  

def calc_eQTL_lm(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    #table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')

    frame = None
    save_path = "/global/scratch/users/ryo10244201/analysis/eQTL_celltype/"

    path_data = "/clusterfs/genomicdata/GTEx/eQTL_files/"

    path_exp = "/global/scratch/users/ryo10244201/analysis/CIBERSORT/data/"

    celltype = input

    tissue = "Whole_Blood"

    try:
        prev = pd.read_csv(save_path + "{}_lm_results.tsv".format(celltype), sep='\t')
        prev_genes = set(prev['gene'])
    except:
        prev_genes = set([])    

    #covariates = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL_covariates/{}.v8.covariates.txt".format(celltype), sep='\t')

    covariates = pd.read_csv("/global/scratch/users/ryo10244201/analysis/Predixcan_cov/{}_covariates.csv".format(tissue), index_col=[0])

    # covariates = covariates.set_index('ID').loc[['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'InferredCov1', 'InferredCov2',
    #              'InferredCov3', 'InferredCov4', 'InferredCov5', 'InferredCov6', 'InferredCov7', 'InferredCov8', 'InferredCov9',
    #               'InferredCov10', 'InferredCov11', 'InferredCov12', 'InferredCov13', 'InferredCov14', 'InferredCov15', 'pcr',
    #               'platform', 'sex']]

    covariates = covariates.set_index('SUBJID')

    covariates = covariates.transpose()

    egenes = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL/{}.v8.egenes.txt.gz".format(tissue), sep='\t', compression='gzip')

    exp_mat = pd.read_csv(path_exp +  "CIBERSORTxHiRes_NA_{}_Window24.txt".format(input), index_col=0, sep='\t')
    def stable_covnert(genesym):
        try:
            return dict_rev[genesym]
        except:
            return None
    exp_mat.index = [stable_covnert(genesym) for genesym in exp_mat.index]

    exp_mat = exp_mat.dropna()
    #exp_mat = exp_mat[exp_mat.columns[3:]]

    genotypes = VariantFile(path_data + "genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz")

    head = str(genotypes.header)
    colums = head[head.rfind('#'):head.rfind("\n")]
    colums = colums.split("\t")

    #print(head)

    lst = []
    lst_all = []

    i = 0
    for gene_id, variant_id, chrom, var_pos in zip(egenes['gene_id'], egenes['variant_id'], egenes['gene_chr'], egenes['variant_pos']):
        print(gene_id)
        if gene_id in prev_genes:
            # sex_temp = prev.loc[prev['gene'] == gene_id]['sex']
            #print(sex_temp)
            continue
        gene_id = gene_id.split(".")[0]
        if not gene_id in exp_mat.index:
            continue
        exp_temp = exp_mat.loc[gene_id]

        exp_temp = exp_temp.dropna()

        if len(exp_temp) == 0:
            continue
        for rec in genotypes.fetch(chrom, var_pos-1, var_pos+1):
            row = str(rec)
            row = row[:-1].split('\t')
            geno_temp = pd.DataFrame(data=[row[9:]], index=['genotype'], columns=colums[9:])
            genotype = row[2]
            #geno_temp = geno_temp.apply(edi_var)
        
        

        exp_frame = pd.DataFrame(data=[exp_temp.values], index=['exp'], columns=exp_temp.index)

        temp_frame = pd.concat([exp_frame, geno_temp, covariates], axis=0, join='inner')

        temp_frame = temp_frame.transpose()

        #print(temp_frame, set(temp_frame['V8']))

        #temp_frame = temp_frame.loc[temp_frame['V8'] == sex]

        #print(temp_frame)

        temp_frame_young = temp_frame.loc[temp_frame['AGE'] < 55]

        temp_frame_old = temp_frame.loc[temp_frame['AGE'] >= 55]

        len_old = temp_frame_old.shape[0]

        len_young = temp_frame_young.shape[0]

        downsample = min(len_old, len_young)

        temp_frame_old = temp_frame_old.sample(downsample)

        temp_frame_young = temp_frame_young.sample(downsample)

        #print(temp_frame_young, temp_frame_old)

        #temp_frame_young, temp_frame_old = split_age_new(temp_frame, table, 55)

        #temp_frame_young = temp_frame_young.transpose()

        temp_frame_young = temp_frame_young.dropna()

        temp_frame_young = temp_frame_young.apply(pd.to_numeric, errors='ignore')

        # mod_young = smf.ols("exp ~ genotype + PC1 + PC2 + PC3 + PC4 + PC5 + InferredCov1 + InferredCov2"
        #  "+ InferredCov3 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9"
        #           "+ InferredCov10 + InferredCov11 + InferredCov12 + InferredCov13 + InferredCov14 + InferredCov15"
        #           "+ pcr + platform + sex", 
        #           data = temp_frame_young)

        temp_frame = temp_frame.dropna()

        temp_frame = temp_frame.apply(pd.to_numeric, errors='ignore')

        # mod_neut = smf.ols("exp ~ genotype + AGE" + "".join([" + V" + str(i) for i in range(1, 19)]), data = temp_frame)

        # res_neut = mod_neut.fit()

        # neut_pval = pd.DataFrame(res_neut.pvalues)
        # neut_pval.index  = "p_val_" + neut_pval.index

        # neut_pval = neut_pval.transpose()

        # neut_slope = pd.DataFrame(res_neut.params)
        # neut_slope.index = "slope_" + neut_slope.index

        # neut_slope = neut_slope.transpose()

        # neut_R2 = pd.DataFrame(data=[res_neut.rsquared], columns=['R2'])

        # neut_meta = pd.DataFrame(data=[[tissue, gene_id, variant_id, 'neut', sex]], columns=["Tissue", "gene", "variant_id","age", 'sex'])

        # neut_frame = pd.concat([neut_meta, neut_slope, neut_pval, neut_R2], axis=1)

        mod_young = smf.ols("exp ~ genotype " + "".join([" + V" + str(i) for i in range(1, 19)]), data = temp_frame_young)

        res_young = mod_young.fit()

        young_pval = pd.DataFrame(res_young.pvalues)
        young_pval.index  = "p_val_" + young_pval.index

        young_pval = young_pval.transpose()

        young_slope = pd.DataFrame(res_young.params)
        young_slope.index = "slope_" + young_slope.index

        young_slope = young_slope.transpose()

        young_R2 = pd.DataFrame(data=[res_young.rsquared], columns=['R2'])

        young_meta = pd.DataFrame(data=[[celltype, gene_id, variant_id, 'young']], columns=["celltype", "gene", "variant_id","age"])

        young_frame = pd.concat([young_meta, young_slope, young_pval, young_R2], axis=1)

        #temp_frame_old = temp_frame_old.transpose()

        temp_frame_old = temp_frame_old.dropna()

        temp_frame_old = temp_frame_old.apply(pd.to_numeric, errors='ignore')

        # mod_old = smf.ols("exp ~ genotype + PC1 + PC2 + PC3 + PC4 + PC5 + InferredCov1 + InferredCov2"
        #  "+ InferredCov3 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9"
        #           "+ InferredCov10 + InferredCov11 + InferredCov12 + InferredCov13 + InferredCov14 + InferredCov15"
        #           "+ pcr + platform + sex", 
        #           data = temp_frame_old)

        mod_old = smf.ols("exp ~ genotype " + "".join([" + V" + str(i) for i in range(1, 19)]), data = temp_frame_old)

        res_old = mod_old.fit()

        old_pval = pd.DataFrame(res_old.pvalues)
        old_pval.index  = "p_val_" + old_pval.index

        old_pval = old_pval.transpose()

        old_slope = pd.DataFrame(res_old.params)
        old_slope.index = "slope_" + old_slope.index

        old_slope = old_slope.transpose()

        old_R2 = pd.DataFrame(data=[res_old.rsquared], columns=['R2'])

        old_meta = pd.DataFrame(data=[[celltype, gene_id, variant_id, 'old']], columns=["celltype", "gene", "variant_id","age"])

        old_frame = pd.concat([old_meta, old_slope, old_pval, old_R2], axis=1)

        frame_temp = pd.concat([young_frame, old_frame
        #, neut_frame
        ], sort=True)

        lst.append(frame_temp)
        lst_all.append(frame_temp)

        if i % 500 == 0:
            print("saving ", i)
            frame_temp = pd.concat(lst)
            if frame is None:
                frame = frame_temp
            else:
                frame = pd.concat([frame_temp, frame])
            lst = []
            frame.to_csv(save_path + "{}_lm_results.tsv".format(celltype), sep='\t')
        i += 1
    if len(lst_all) > len(prev_genes):
        frame = pd.concat(lst_all)
        frame.to_csv(save_path + "{}_lm_results.tsv".format(celltype), sep='\t')

def edi_var(var, id):
    zero = id[0]
    one = id[1]
    if var == "0|0":
        return 0
    elif var == "0|1":
        return 1
    elif var == "1|0":
        return 1
    else:
        return 2

if args.test == True:
    print("testing")
    celltypes = [
        "Bcells",
        "Monocytes",
        "NKcells",
        "NKTcells",
        "TcellsCD4",
        "TcellsCD8"
        ]
    for celltype in celltypes:
        calc_eQTL_lm(celltype ,"")

celltype = snakemake.input[0].split("_")[-2]

calc_eQTL_lm(celltype, snakemake.output[0])
