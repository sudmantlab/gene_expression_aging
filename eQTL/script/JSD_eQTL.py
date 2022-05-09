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

np.random.seed(0)

parser = argparse.ArgumentParser(description='Snakemake script for calculating adjacent gene Spearman correlation.')
parser.add_argument('-t', '--test', action='store_true', 
                    default=False,
                    help='testing script with single Blood input')

parser.add_argument('-f', '--frame', action='store_true', 
                    default=False,
                    help='saving frame with single Blood input')

args = parser.parse_args()

def calc_eQTL_lm(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')

    frame = None
    save_path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/analysis/JSD_eQTL_new/"

    path_data = "/clusterfs/genomicdata/GTEx/eQTL_files/"

    tissue = input

    try:
        prev = pd.read_csv(save_path + "{}_lm_results_new.tsv".format(tissue), sep='\t')
        prev_genes = set(prev['gene'])
    except:
        prev_genes = set([])    

    #covariates = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL_covariates/{}.v8.covariates.txt".format(tissue), sep='\t')

    covariates = pd.read_csv("/global/scratch/users/ryo10244201/analysis/Predixcan_cov/{}_covariates.csv".format(tissue), index_col=[0])

    # covariates = covariates.set_index('ID').loc[['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'InferredCov1', 'InferredCov2',
    #              'InferredCov3', 'InferredCov4', 'InferredCov5', 'InferredCov6', 'InferredCov7', 'InferredCov8', 'InferredCov9',
    #               'InferredCov10', 'InferredCov11', 'InferredCov12', 'InferredCov13', 'InferredCov14', 'InferredCov15', 'pcr',
    #               'platform', 'sex']]

    covariates = covariates.set_index('SUBJID')

    covariates = covariates.transpose()

    egenes = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL/{}.v8.egenes.txt.gz".format(tissue), sep='\t', compression='gzip')

    exp_mat = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL_expression_matrices/{}.v8.normalized_expression.bed.gz".format(tissue), sep='\t', compression='gzip')
    exp_mat = exp_mat.set_index("gene_id")
    exp_mat = exp_mat[exp_mat.columns[3:]]

    genotypes = VariantFile(path_data + "genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz")

    head = str(genotypes.header)
    colums = head[head.rfind('#'):head.rfind("\n")]
    colums = colums.split("\t")

    lst = []

    i = 0
    for gene_id, variant_id, chrom, var_pos in zip(egenes['gene_id'], egenes['variant_id'], egenes['gene_chr'], egenes['variant_pos']):
        if gene_id in prev_genes:
            continue
        for rec in genotypes.fetch(chrom, var_pos-1, var_pos+1):
            row = str(rec)
            row = row[:-1].split('\t')
            geno_temp = pd.DataFrame(data=[row[9:]], index=['genotype'], columns=colums[9:])
            genotype = row[2]
        
        exp_temp = exp_mat.loc[gene_id]

        exp_frame = pd.DataFrame(data=[exp_temp.values], index=['exp'], columns=exp_temp.index)

        temp_frame = pd.concat([exp_frame, geno_temp, covariates], axis=0, join='inner')

        temp_frame = temp_frame.transpose()

        temp_frame_young = temp_frame.loc[temp_frame['AGE'] < 55]

        temp_frame_old = temp_frame.loc[temp_frame['AGE'] >= 55]

        len_old = temp_frame_old.shape[0]

        len_young = temp_frame_young.shape[0]

        downsample = min(len_old, len_young)

        temp_frame_old = temp_frame_old.sample(downsample)

        temp_frame_young = temp_frame_young.sample(downsample)

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

        mod_neut = smf.ols("exp ~ genotype + AGE" + "".join([" + V" + str(i) for i in range(1, 19)]), data = temp_frame)

        res_neut = mod_neut.fit()

        neut_pval = pd.DataFrame(res_neut.pvalues)
        neut_pval.index  = "p_val_" + neut_pval.index

        neut_pval = neut_pval.transpose()

        neut_slope = pd.DataFrame(res_neut.params)
        neut_slope.index = "slope_" + neut_slope.index

        neut_slope = neut_slope.transpose()

        neut_R2 = pd.DataFrame(data=[res_neut.rsquared], columns=['R2'])

        neut_meta = pd.DataFrame(data=[[tissue, gene_id, variant_id, 'neut']], columns=["Tissue", "gene", "variant_id","age"])

        neut_frame = pd.concat([neut_meta, neut_slope, neut_pval, neut_R2], axis=1)

        mod_young = smf.ols("exp ~ genotype " + "".join([" + V" + str(i) for i in range(1, 19)]), data = temp_frame_young)

        res_young = mod_young.fit()

        young_pval = pd.DataFrame(res_young.pvalues)
        young_pval.index  = "p_val_" + young_pval.index

        young_pval = young_pval.transpose()

        young_slope = pd.DataFrame(res_young.params)
        young_slope.index = "slope_" + young_slope.index

        young_slope = young_slope.transpose()

        young_R2 = pd.DataFrame(data=[res_young.rsquared], columns=['R2'])

        young_meta = pd.DataFrame(data=[[tissue, gene_id, variant_id, 'young']], columns=["Tissue", "gene", "variant_id","age"])

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

        old_meta = pd.DataFrame(data=[[tissue, gene_id, variant_id, 'old']], columns=["Tissue", "gene", "variant_id","age"])

        old_frame = pd.concat([old_meta, old_slope, old_pval, old_R2], axis=1)

        frame_temp = pd.concat([young_frame, old_frame, neut_frame])

        lst.append(frame_temp)

        if i % 5000 == 0:
            if frame is None:
                frame = frame_temp
            else:
                frame = pd.concat([frame_temp, frame])
            frame.to_csv(save_path + "{}_lm_results_new.tsv".format(tissue), sep='\t')
        i += 1
    if len(lst) > len(prev_genes):
        frame = pd.concat(lst)
        frame.to_csv(save_path + "{}_lm_results_new.tsv".format(tissue), sep='\t')

def calc_eQTL_lm_save(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    #table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')

    frame = None
    save_path = "/global/scratch/users/ryo10244201/analysis/JSD_eQTL/"

    path_data = "/clusterfs/genomicdata/GTEx/eQTL_files/"

    tissue = input

    egenes = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL/{}.v8.egenes.txt.gz".format(tissue), sep='\t', compression='gzip')

    exp_mat = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL_expression_matrices/{}.v8.normalized_expression.bed.gz".format(tissue), sep='\t', compression='gzip')
    exp_mat = exp_mat.set_index("gene_id")
    exp_mat = exp_mat[exp_mat.columns[3:]]

    genotypes = VariantFile(path_data + "genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz")

    head = str(genotypes.header)
    colums = head[head.rfind('#'):head.rfind("\n")]
    colums = colums.split("\t")

    lst = []

    i = 0
    for gene_id, variant_id, chrom, var_pos in zip(egenes['gene_id'], egenes['variant_id'], egenes['gene_chr'], egenes['variant_pos']):
        for rec in genotypes.fetch(chrom, var_pos-1, var_pos+1):
            row = str(rec)
            row = row[:-1].split('\t')
            geno_temp = pd.DataFrame(data=[row[9:]], index=[gene_id + '_geno'], columns=colums[9:])
            genotype = row[2]
        
        print(geno_temp)

        ids = variant_id.split("_")

        var_id = (ids[2], ids[3])

        geno_temp = geno_temp.applymap(lambda x : edi_var(x, var_id))
        print(geno_temp)

        exp_temp = exp_mat.loc[gene_id]

        exp_frame = pd.DataFrame(data=[exp_temp.values], index=[gene_id + '_exp'], columns=exp_temp.index)

        temp_frame = pd.concat([exp_frame, geno_temp])


        #temp_frame['group'] = [0, gene_id]

        lst.append(temp_frame)



        # if frame is None:
        #     frame = temp_frame
        # else:
        #     frame = pd.concat([temp_frame, frame])
        i += 1
    frame = pd.concat(lst)
    frame = frame.transpose()
    frame.to_csv(save_path + "{}_lm_dist_4.tsv".format(tissue), sep='\t')

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
    tissues = [ "Esophagus_Gastroesophageal_Junction"]
    for tissue in tissues:
        calc_eQTL_lm(tissue ,"")

if args.frame == True:
    print("saving frame")
    calc_eQTL_lm_save("Heart_Atrial_Appendage", "")

calc_eQTL_lm(snakemake.input[0][snakemake.input[0].rfind('/') + 1:snakemake.input[0].find('.')], snakemake.output[0])
