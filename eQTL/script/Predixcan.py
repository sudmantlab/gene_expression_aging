import numpy as np
import pandas as pd
import os

path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/analysis/Predixcan_eQTL/"

path_data = "/clusterfs/genomicdata/GTEx/eQTL_files/"

tissue = 'Whole_Blood'

data = np.load(path + "r2_{}_old.npy".format(tissue))

tissues_b9 = ["Lung",
    "Artery_Tibial",
    "Heart_Left_Ventricle",
    "Nerve_Tibial",
    "Thyroid",
    "Adipose_Subcutaneous",
    "Skin_Sun_Exposed_Lower_leg",
    "Muscle_Skeletal",
    "Whole_Blood"]

tissues_b27 = ['Whole_Blood','Stomach',
         'Colon_Sigmoid',
         'Esophagus_Gastroesophageal_Junction',
         'Colon_Transverse',
         'Artery_Aorta',
         'Heart_Atrial_Appendage',
         'Breast_Mammary_Tissue',
         'Prostate',
         'Heart_Left_Ventricle',
         'Esophagus_Muscularis',
         'Adipose_Visceral_Omentum',
         'Esophagus_Mucosa',
         'Lung',
         'Skin_Not_Sun_Exposed_Suprapubic',
         'Nerve_Tibial',
         'Thyroid',
         'Adipose_Subcutaneous',
         'Artery_Tibial',
         'Testis',
         'Skin_Sun_Exposed_Lower_leg',
         'Muscle_Skeletal',
         'Pancreas',
         'Liver',
         'Artery_Coronary',
         'Adrenal_Gland',
         'Cells_Cultured_fibroblasts']

frame = pd.DataFrame(columns=['Tissue', 'age', 'gene', 'R2'], data=[])

for tissue in tissues_b9:
    exp_mat = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL_expression_matrices/{}.v8.normalized_expression.bed.gz".format(tissue), sep='\t', compression='gzip')
    exp_mat = exp_mat.set_index("gene_id")
    print(exp_mat.index)
    data_old = np.load(path + "r2_{}_old.npy".format(tissue))
    sz = len(data_old)
    #break
    frame_old = pd.DataFrame(columns=['Tissue', 'age', 'gene', 'R2'], data=np.array([[tissue]*sz, ['old']*sz, exp_mat.index.values, data_old]).T )
    data_young = np.load(path + "r2_{}_young.npy".format(tissue))
    sz = len(data_young)
    #break
    frame_young = pd.DataFrame(columns=['Tissue', 'age', 'gene', 'R2'], data=np.array([[tissue]*sz, ['young']*sz, exp_mat.index.values, data_young]).T )

    frame = pd.concat([frame, frame_old, frame_young], axis=0)

frame.to_csv(path + "r2_b9.tsv", sep='\t')

for tissue in tissues_b27:
    exp_mat = pd.read_csv(path_data + "GTEx_Analysis_v8_eQTL_expression_matrices/{}.v8.normalized_expression.bed.gz".format(tissue), sep='\t', compression='gzip')
    exp_mat = exp_mat.set_index("gene_id")
    print(exp_mat.index)
    data_old = np.load(path + "r2_{}_old.npy".format(tissue))
    sz = len(data_old)
    #break
    frame_old = pd.DataFrame(columns=['Tissue', 'age', 'gene', 'R2'], data=np.array([[tissue]*sz, ['old']*sz, exp_mat.index.values, data_old]).T )
    data_young = np.load(path + "r2_{}_young.npy".format(tissue))
    sz = len(data_young)
    #break
    frame_young = pd.DataFrame(columns=['Tissue', 'age', 'gene', 'R2'], data=np.array([[tissue]*sz, ['young']*sz, exp_mat.index.values, data_young]).T )

    frame = pd.concat([frame, frame_old, frame_young], axis=0)

frame.to_csv(path + "r2_b27.tsv", sep='\t')
