import os
import pandas as pd
from scipy.stats import entropy, ttest_ind_from_stats
import numpy as np
from scipy.spatial import distance as dist
import argparse
from JSD_utils import *

parser = argparse.ArgumentParser(description='Script for calculating variance distribution.')
parser.add_argument('-t', '--test', action='store_true', 
                    default=False,
                    help='testing script with single input')

args = parser.parse_args()

def calc_binned_JSD(input, output):
    path = "/global/home/users/shenghuanjie/projects/ge_variance/datasets/gtex_new/input/input_by_tissue/"
    table = pd.read_csv(path + "../phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt.gz", sep='\t', compression='gzip')
    frame = pd.DataFrame(columns=["Tissue", "metric", "bin", "mean"])
    dist_frame = pd.DataFrame(columns=["Tissue", "metric", "bin", "age", "value"])
    save_path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/analysis/JSD_binned/"
    temp = pd.read_csv(input, sep="\t")

    tissue = input[input.rfind('/') + 1:][14:-4]

    binned_frames = split_age_bins(temp, table, 6)

    #print(temp)
    #Jensen-Shanon Divergence
    print("JSD start")

    i = 1
    for binned_frame in binned_frames:
        if binned_frame.size == 0:

            print("no sample in this bin")
            continue
        try:
            print(binned_frame)
            jsd_mean, jsd_dist, age = pairwise_JSD_age(binned_frame, table)
        except ValueError:
            print("no sample in this bin")
            continue
        frame = frame.append(pd.DataFrame(data=[[tissue, "JSD" ,i, jsd_mean]], columns=frame.columns))
        print(dist_frame)
        dist_frame = pd.concat([dist_frame, pd.DataFrame(data=np.array([np.repeat(tissue, len(jsd_dist)), np.repeat("JSD", len(jsd_dist)),\
             np.repeat(i, len(jsd_dist)), age, jsd_dist]).T, columns=dist_frame.columns)], axis=0)
        i += 1

    frame.to_csv(save_path + tissue + "_summary.csv")
    dist_frame.to_csv(output)

if args.test:
    print("testing")
    calc_binned_JSD("/global/scratch/shenghuanjie/datasets/gtex_new/input/input_by_tissue/GTEx_Analysis-Cells_Culturedfibroblasts.tsv",\
        "/global/scratch/ryo10244201/analysis/JSD_binned/Cells_Culturedfibroblasts_dist.csv")

print(snakemake.input[0], snakemake.output[0])
calc_binned_JSD(snakemake.input[0], snakemake.output[0])
