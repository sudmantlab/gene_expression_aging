import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues
import time

def app_combine_pvalues(pvalues):
    np_pvalues = np.array(pvalues)
    np_pvalues = np_pvalues[~np.isnan(pvalues)]
    return combine_pvalues(np_pvalues)[1]
start = time.time()
print(time.time())
save_path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/analysis/JSD_eQTL_new/"

tissues = ["Whole_Blood",
        "Stomach",
        "Colon_Sigmoid",
        "Esophagus_Gastroesophageal_Junction",
        "Colon_Transverse",
        "Artery_Aorta",
        "Heart_Atrial_Appendage",
        "Breast_Mammary_Tissue",
        "Prostate",
        "Heart_Left_Ventricle",
         "Esophagus_Muscularis",
         "Adipose_Visceral_Omentum",
         "Esophagus_Mucosa",
         "Lung",
         "Skin_Not_Sun_Exposed_Suprapubic",
         "Nerve_Tibial",
         "Thyroid",
         "Adipose_Subcutaneous",
         "Artery_Tibial",
         "Testis",
         "Skin_Sun_Exposed_Lower_leg",
         "Muscle_Skeletal",
         "Pancreas",
         "Liver",
         "Artery_Coronary",
         "Adrenal_Gland",
         "Cells_Cultured_fibroblasts"]

tissue_path = [save_path + tissue + '_lm_results_new.tsv' for tissue in tissues]

temp = pd.concat([pd.read_csv(f, sep='\t') for f in tissue_path])[["R2","Tissue", "gene", "age", "p_val_AGE", "slope_AGE","p_val_genotype[T.0|1]", "p_val_genotype[T.1|0]", "p_val_genotype[T.1|1]"]]

temp['pval_combined'] = temp[["p_val_genotype[T.0|1]", "p_val_genotype[T.1|0]", "p_val_genotype[T.1|1]"]].apply(app_combine_pvalues, axis=1)

print(temp['pval_combined'])

temp.to_csv(save_path + "All_Tissue_lm_results_new.tsv", sep='\t')

print(time.time() - start)

# frame_all = pd.read_csv(save_path + "All_lm_results_2.tsv", sep='\t')

# frame_all = pd.concat([temp, frame_all])

# frame_all.to_csv(save_path + "All_lm_results_2.tsv", sep='\t')


