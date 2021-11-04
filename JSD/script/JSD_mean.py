import os
import pandas as pd
from scipy.stats import entropy, ttest_ind_from_stats
import numpy as np
from scipy.spatial import distance as dist
from JSD_utils import *

def calc_JSD_mean(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')
    frame = pd.DataFrame(columns=["Tissue", "gene" , "young", "old"])
    save_path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/analysis/JSD_All/"
    temp = pd.read_csv(input, sep="\t")

    print(temp)

    temp = temp.sort_values()

    tissue = input[input.rfind('/') + 1:][14:-4]

    young, old= split_age_new(temp, table, 55)

    mean_young = np.mean(young)

    jsd_old, jsd_old_dist  = pairwise_JSD(old)

    frame = frame.append(pd.DataFrame(data=[[tissue, mean_young, mean_old]],
    columns=frame.columns))

    frame.to_csv(output)

calc_JSD_mean("/global/scratch/shenghuanjie/datasets/gtex_new/input/input_by_tissue/GTEx_Analysis-Esophagus_Mucosa.tsv",\
        "/global/scratch/ryo10244201/analysis/JSD_All/Esophagus_Mucosa_JSD_result.csv")
print(snakemake.input[0], snakemake.output[0])
calc_JSD_mean(snakemake.input[0], snakemake.output[0])
