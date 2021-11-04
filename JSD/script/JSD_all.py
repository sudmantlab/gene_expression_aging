import os
import pandas as pd
from scipy.stats import entropy, ttest_ind_from_stats
import numpy as np
from scipy.spatial import distance as dist
from JSD_utils import *

def calc_JSD(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')
    frame = pd.DataFrame(columns=["Tissue", "p-value" , "young", "old"])
    save_path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/analysis/JSD_All/"
    temp = pd.read_csv(input, sep="\t")

    tissue = input[input.rfind('/') + 1:][14:-4]

    young, old= split_age_new(temp, table, 55)

    jsd_young, jsd_young_dist = pairwise_JSD(young)

    jsd_old, jsd_old_dist  = pairwise_JSD(old)

    p = ttest_ind_from_stats(jsd_old, np.std(jsd_old_dist), jsd_old_dist.size, 
    jsd_young, np.std(jsd_young_dist), jsd_young_dist.size, False)
    frame = frame.append(pd.DataFrame(data=[[tissue, p[1], jsd_young, jsd_old]],
    columns=frame.columns))

    frame.to_csv(output)
    np.savetxt(save_path + tissue + "_young_dist.csv", jsd_young_dist, delimiter=',')
    np.savetxt(save_path + tissue + "_old_dist.csv", jsd_old_dist, delimiter=',')

calc_JSD("/global/scratch/shenghuanjie/datasets/gtex_new/input/input_by_tissue/GTEx_Analysis-Esophagus_Mucosa.tsv",\
        "/global/scratch/ryo10244201/analysis/JSD_All/Esophagus_Mucosa_JSD_result.csv")
print(snakemake.input[0], snakemake.output[0])
calc_JSD(snakemake.input[0], snakemake.output[0])
