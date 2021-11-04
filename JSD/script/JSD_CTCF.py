import os
import pandas as pd
from scipy.stats import entropy, ttest_ind_from_stats, pearsonr, spearmanr
import numpy as np
from scipy.spatial import distance as dist
import time
from JSD_utils import age_search, pairwise_JSD
import argparse

parser = argparse.ArgumentParser(description='Snakemake script for calculating JSD between / not between CTCF region.')
parser.add_argument('-t', '--test', action='store_true', 
                    default=False,
                    help='testing script with single input')
args = parser.parse_args()

def gene_mapping_table(table, df):
    dic = {}
    somelst = [str(num) for num in range(1, 23)]
    somelst.extend(['X'])
    for sex in somelst[::-1]:
        c = [val for val in  table.loc[table['chromosome_name'] == sex]['Name#Description']]
        b = [val in c for val in df["Name#Description"].values]
        a = df.loc[b]
        if len(a.index) == 0:
            continue
        dic[sex] = a
    return dic

def split_age(df, table, threshold):
    """Return split of temp dataframe given table of sample info according to threshold"""
    agelist = age_search(df.columns.values[1:], table)
    row = df["Name#Description"]
    df = df.drop(df.columns[0], axis=1)
    df = df.transpose()
    young = df.loc[agelist < threshold]
    old = df.loc[agelist >= threshold]
    return young.transpose(), old.transpose(), row

def geneset_contains_CTCF(anno, gene_table):
    somelst = [str(num) for num in range(1, 23)]
    somelst.extend(['X'])
    ret = {}
    for chrom in somelst:
        tic = time.time()
        print("start anno subset " + str(tic))
        anno_temp = anno.loc[anno['chrom'] == 'chr' + chrom]
        g_t_temp = gene_table.loc[gene_table['chromosome_name'] == chrom][['Name#Description', 'start_position']]
        print("anno subset: " + str(time.time() - tic))
        i, j = 0, 1
        lst = []
        while i < len(anno_temp) or j < len(anno_temp):
            region1 = anno.iloc[i]
            region2 = anno.iloc[j]
            window_size = region2['start'] - region1['end'] // 2
            i += 1
            j += 1
            if window_size < 500:
                print('skipping because of window size')
                continue
            between = g_t_temp.loc[(g_t_temp['start_position'] > region1['end']) &
             (g_t_temp['start_position'] < region2['start'])]['Name#Description']
            center_reg1 = g_t_temp.loc[(g_t_temp['start_position'] > region1['start'] - window_size )&
              (g_t_temp['start_position'] < region1['end'] + window_size)]['Name#Description']
            center_reg2 = g_t_temp.loc[(g_t_temp['start_position'] > region2['start'] - window_size) &
              (g_t_temp['start_position'] < region2['end'] + window_size)]['Name#Description']
            if len(between) < 3 or len(center_reg1) < 3 or len(center_reg2) < 3:
                #print('skipping because there is no gene')
                continue
            lst.append((between.values, center_reg1.values, center_reg2.values))
            if i % 100 == 0:
                print(i, j, "of", len(anno_temp))
                #print(lst)
        tic = time.time()
        print("start anno geneset " + str(tic))
        ret[chrom] = lst
        print("anno geneset: " + str(time.time() - tic))
    return ret

def calc_JSD_CTCF(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')

    path_anno = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/seg_anno.tsv"
    anno = pd.read_csv(path_anno, sep='\t')

    path_gene = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/"
    gene_table = pd.read_csv(path_gene + 'psh-genelist_ensembl.tsv', sep='\t')

    frame = pd.DataFrame(columns=["Tissue","Chromosome", "contains_CTCF", "young", "old", "p_value"])

    save_path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/JSD_CTCF/"
    temps = pd.read_csv(input, sep="\t")

    i = 0

    gene_table.sort_values(by=['chromosome_name', 'start_position'], inplace=True)
    dic = gene_mapping_table(gene_table, temps)
    tissue = input[input.rfind('/') + 1:][14:-4]
    
    tic = time.time()
    print("start anno geneset " + str(tic))
    anno_dic = geneset_contains_CTCF(anno, gene_table)
    print("anno geneset: " + str(time.time() - tic))
    for chrom in dic:
        
        young, old, row = split_age(dic[chrom], table, 55)
        
        for bet, rg1, rg2 in anno_dic[chrom]:
            bet, rg1, rg2 = [gene in bet for gene in dic[chrom]['Name#Description']], [gene in rg1 for gene in dic[chrom]['Name#Description']], [gene in rg2 for gene in dic[chrom]['Name#Description']]
            young_bet, young_rg1, young_rg2 = young.loc[bet], young.loc[rg1], young.loc[rg2]
            old_bet, old_rg1, old_rg2 = old.loc[bet], old.loc[rg1], old.loc[rg2]
            print("{} out of {}".format(i, len(anno_dic[chrom])))
            if len(young_bet) == 0 or len(young_rg1) == 0 or len(young_rg2) == 0 or len(old_bet) == 0 or len(old_rg1) == 0 or len(old_rg2) == 0:
                print('skipping because there is no gene')
                continue
            young_bet_JSD, young_bet_dist = pairwise_JSD(young_bet)
            old_bet_JSD, old_bet_dist = pairwise_JSD(old_bet)
            p_bet = ttest_ind_from_stats(old_bet_JSD, np.std(old_bet_dist), old_bet_dist.size, 
            young_bet_JSD, np.std(young_bet_dist), young_bet_dist.size, False)
            frame = frame.append({"Tissue": tissue, "Chromosome": chrom, "contains_CTCF": False,
             "young": young_bet_JSD, "old": old_bet_JSD, "p_value": p_bet[1]}, ignore_index=True)

            young_rg1_JSD, young_rg1_dist = pairwise_JSD(young_rg1)
            old_rg1_JSD, old_rg1_dist = pairwise_JSD(old_rg1)
            p_rg1 = ttest_ind_from_stats(old_rg1_JSD, np.std(old_rg1_dist), old_rg1_dist.size, 
            young_rg1_JSD, np.std(young_rg1_dist), young_rg1_dist.size, False)
            frame = frame.append({"Tissue": tissue, "Chromosome": chrom, "contains_CTCF": True,
             "young": young_rg1_JSD, "old": old_rg1_JSD, "p_value": p_rg1[1]}, ignore_index=True)
            
            young_rg2_JSD, young_rg2_dist = pairwise_JSD(young_rg2)
            old_rg2_JSD, old_rg2_dist = pairwise_JSD(old_rg2)
            p_rg2 = ttest_ind_from_stats(old_rg2_JSD, np.std(old_rg2_dist), old_rg2_dist.size, 
            young_rg2_JSD, np.std(young_rg2_dist), young_rg2_dist.size, False)
            frame = frame.append({"Tissue": tissue, "Chromosome": chrom, "contains_CTCF": True,
             "young": young_rg2_JSD, "old": old_rg2_JSD, "p_value": p_rg2[1]}, ignore_index=True)
            if i % 10 == 0:
                print('save at 10')
                frame.to_csv(save_path + input[input.rfind('/') + 1:][14:-4] + "_JSD_CTCF_result.tsv", sep='\t')
            i += 1
    print('result frame', frame)
    frame.to_csv(save_path + input[input.rfind('/') + 1:][14:-4] + "_JSD_CTCF_result.tsv", sep='\t')

if args.test == True:
    print("testing")
    calc_JSD_CTCF("/global/scratch/shenghuanjie/datasets/gtex_new/input/input_by_tissue/GTEx_Analysis-WholeBlood.tsv" ,"")


calc_JSD_CTCF(snakemake.input[0], snakemake.output[0])

