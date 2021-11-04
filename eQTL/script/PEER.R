library(peer)
library(vroom)
library(dplyr)
library(effectsize)

path = "/clusterfs/genomicdata/GTEx/eQTL_files/GTEx_Analysis_v8_eQTL_expression_matrices/"

path_sc = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/data/misc/"



gt50 = c(
    # 'Whole_Blood',
    # 'Stomach',
    #      'Colon_Sigmoid',
    #      'Esophagus_Gastroesophageal_Junction',
    #      'Colon_Transverse',
    #      'Artery_Aorta',
    #      'Heart_Atrial_Appendage',
        #  'Breast_Mammary_Tissue',
        #  'Prostate',
        #  'Heart_Left_Ventricle',
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
         'Cells_Cultured_fibroblasts'
         )

for (tis in gt50) {
    frame = vroom(paste0(path_sc, '/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt') )

    age = frame %>% select(SUBJID, AGE) %>% unique() %>% select(AGE, SUBJID)

    age$AGE <- as.numeric(age$AGE)

    fr_expr = vroom(paste0(path, tis, ".v8.normalized_expression.bed.gz"))

    fr_expr_t = as.data.frame(t(as.matrix( fr_expr %>% select(-c(`#chr`, "start", "end", "gene_id"))
    )))

    names(fr_expr_t) = fr_expr$gene_id

    fr_cov = vroom(paste0(path, "../GTEx_Analysis_v8_eQTL_covariates/",tis, '.v8.covariates.txt') )

    ID_col = fr_cov['ID']

    new_fr_cov = subset(fr_cov, select = -ID)

    rownames(new_fr_cov) = ID_col$ID

    fr_cov_t = as.data.frame(t(as.matrix(new_fr_cov)))

    fr_cov_t$SUBJID = rownames(fr_cov_t)

    new_merged = merge(fr_cov_t, age, by = 'SUBJID')

    new_merged_cov = new_merged %>% select(PC1, PC2, PC3, PC4, PC5, pcr, platform, sex)

    fr_expr_t$SUBJID = fr_cov_t$SUBJID

    expr_age = merge(fr_expr_t, age, by='SUBJID')

    adjust_expr = adjust(expr_age %>% select(-c('SUBJID')), effect = 'AGE', exclude= 'AGE')

    adjust_expr_save = adjust_expr

    adjust_expr_save$SUBJID = expr_age$SUBJID

    write.csv(adjust_expr_save, paste0(path_sc, "../../analysis/Predixcan_cov/", tis, "_adjust_exp.csv"))

    adjust_expr = adjust_expr %>% select(-c('AGE'))

    model = PEER()

    PEER_setPhenoMean(model, as.matrix(adjust_expr))

    PEER_setNk(model,10)

    PEER_setCovariates(model, as.matrix(new_merged_cov))

    PEER_update(model)

    factors = PEER_getX(model)

    cof_frame = as.data.frame(factors)

    cof_frame$AGE = new_merged$AGE

    cof_frame$SUBJID = new_merged$SUBJID

    write.csv(cof_frame, paste0(path_sc, "../../analysis/Predixcan_cov/", tis, "_covariates.csv"))


    line = lm(formula = AGE~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18,
    data = cof_frame)

}

fr_expr = vroom(paste0(path, "Whole_Blood.v8.normalized_expression.bed.gz"))

fr_expr_t = as.data.frame(t(as.matrix( fr_expr %>% select(-c(`#chr`, "start", "end", "gene_id"))
)))

names(fr_expr_t) = fr_expr$gene_id

fr_cov = vroom(paste0(path_sc, 'Whole_Blood.v8.covariates.txt') )

ID_col = fr_cov['ID']

new_fr_cov = subset(fr_cov, select = -ID)

rownames(new_fr_cov) = ID_col$ID

fr_cov_t = as.data.frame(t(as.matrix(new_fr_cov)))

fr_cov_t$SUBJID = rownames(fr_cov_t)

new_merged = merge(fr_cov_t, age, by = 'SUBJID')

new_merged_cov = new_merged %>% select(PC1, PC2, PC3, PC4, PC5, pcr, platform, sex)

fr_expr_t$SUBJID = fr_cov_t$SUBJID

expr_age = merge(fr_expr_t, age, by='SUBJID')

adjust_expr = adjust(expr_age %>% select(-c('SUBJID')), effect = 'AGE', exclude= 'AGE')

adjust_expr = adjust_expr %>% select(-c('AGE'))

model = PEER()

PEER_setPhenoMean(model, as.matrix(adjust_expr))

PEER_setNk(model,10)

PEER_setCovariates(model, as.matrix(new_merged_cov))

PEER_update(model)

factors = PEER_getX(model)

cof_frame = as.data.frame(factors)

cof_frame$AGE = new_merged$AGE


line = lm(formula = AGE~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18,
 data = cof_frame)


