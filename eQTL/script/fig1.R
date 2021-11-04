library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
# library(tidyr)
library(scales)
library(ggrepel)
library(ggsci)
library(vroom)
library(dlookr)
library(ggpubr)
library(purrr)
library(broom)
# setwd("Y:/")
library(grid)
library(relaimpo)
setwd("Z:/sudmant_lab/gene_expression_variance/Rscript")
source('../Rscript/cleanSubType.R')
library(rstan)
library(BisRNA)
library(latex2exp)
library(PNWColors)
library(stringr)
library(data.table)

source('../Rscript/cleanSubType.R')

path_temp <- 'G:/Shared drives/sudmantlab/projects/agingGeneRegulation/figures/figure1/fig1_'

samp_table <- vroom("Y:/data/misc/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Info.txt")

samp_table = samp_table %>%
  dplyr::select(SUBJID, AGE) %>%
  unique() 

data_path <- "Y:/analysis/JSD_eQTL/Whole_Blood_lm_dist_4.tsv"

frame <- vroom(data_path)

genes = c("ENSG00000137106.17", "ENSG00000258289.8",
          "ENSG00000273295.1",  "ENSG00000214425.7", 
          "ENSG00000178952.10", "ENSG00000074803.17",
          "ENSG00000108342.12")

frame %>% dplyr::select(starts_with('ENSG00000108342'))

panel = function(gene) {
  frame_gene = frame %>%
    dplyr::select("...1", starts_with(gene)) %>%
    dplyr::mutate(SUBJID = ...1, exp = !!as.name(paste0(gene, '_exp')), genotype = !!as.name(paste0(gene, '_geno')))  %>%
    filter(!is.na(exp)) %>%
    mutate(genotype = str_replace(genotype, '_', '/'))
    group_by(genotype)
  
  test_t = inner_join(frame_gene, samp_table, by='SUBJID')
  
  g1 = ggplot(test_t %>%
                filter(!is.na(exp)) %>%
                mutate(age = ifelse(AGE > 55, 'old', 'young'))
  ) +
    stat_summary(aes(y = exp, x = fct_rev(age), col= age), size = 1, alpha = 0.8) +
    geom_jitter(aes(y = exp, x = fct_rev(age), col= age) , alpha = 0.3) +
    facet_wrap(~genotype, nrow = 1) + 
    labs(
      x=TeX('genotype'),
      y=TeX('expression')
    )+
    theme_bw(base_size=12)+
    theme(
      legend.position = "None",
      legend.box = "vertical",
      text = element_text(size=10),
      axis.title = element_text(size=12),
      axis.text = element_text(size=12),
      strip.background = element_blank(),
      strip.text = element_text(size=12)
    ) + 
    scale_color_brewer(
      palette = 'Set1'
    )
  return(g1)
}

genes_plot = c("ENSG00000108342.12", "ENSG00000178952.10", "ENSG00000258289.8", "ENSG00000214425.7")

g1.1 = panel("ENSG00000108342.12")

g1.1

g1.2 = panel("ENSG00000178952.10")

g1.3 = panel("ENSG00000258289.8")

g1.4 = panel("ENSG00000214425.7")

gene_ids = c("CSF3", "TUFM", "CHURC1", "LRRC37A4P")

g1 = ggarrange(g1.1, g1.2, g1.3, g1.4, ncol = 2, nrow=2,
               label.y = 1.02,
               font.label = list(size = 12, color = "black",face="bold.italic", family = NULL), 
               labels = gene_ids)

g1

ggsave(paste0(path_temp,'gene_examples.png'), plot = g1, width = 14, height = 12, scale=0.5, dpi=600)

str_split(names(frame), '_', 1)

gene_all = str_extract(names(frame), "(.*)(?=_)") %>% 
  unique()

gene_all = na.omit(gene_all)

full_frame = data.frame(genotype=character(0),
                        age=character(0),
                        mu_g=numeric(0),
                        gene=character(0)
)

for (gene in gene_all) {
  frame_gene = frame %>%
    dplyr::select("...1", starts_with(gene)) %>%
    dplyr::mutate(SUBJID = ...1, exp = !!as.name(paste0(gene, '_exp')), genotype = !!as.name(paste0(gene, '_geno')))  %>%
    filter(!is.na(exp)) %>%
    group_by(genotype)
  
  test_t = inner_join(frame_gene, samp_table, by='SUBJID')
  
  s = test_t %>%
    mutate(age = ifelse(AGE > 55, 'old', 'young')) %>%
    group_by(genotype, age) %>%
    dplyr::summarise(mu_g = mean(exp))
  s$gene = gene
  
  full_frame = rbind(full_frame, s)
}

write.csv2(full_frame, "../analysis/meam_genotype_blood.csv")

full_frame_new = full_frame %>%
  pivot_wider(names_from = age, values_from=mu_g) %>%
  mutate(delta = young - old) %>%
  mutate('trend' = ifelse(delta > 0, 'young_up', 'old_up')) %>%
  filter(!genotype %in% c('A_A', 'G_G', 'C_C', 'T_T'))

full_frame_new %>%
  group_by(trend) %>%
  dplyr::summarise(n=n())

ggplot(full_frame %>% filter(!genotype %in% c('A_A', 'G_G', 'C_C', 'T_T'))) + geom_boxplot(aes(x = age, y=mu_g))

mean(na.omit(full_frame_new$delta))

t.test(full_frame_new$old, full_frame_new$young, paired = TRUE, alternative = "two.sided")

t.test(na.omit(full_frame_new$delta))

frame_gene = frame %>%
  dplyr::select("...1", starts_with(gene)) %>%
  dplyr::mutate(SUBJID = ...1, exp = !!as.name(paste0(gene, '_exp')), genotype = !!as.name(paste0(gene, '_geno')))  %>%
  filter(!is.na(exp)) %>%
  group_by(genotype)



test_t = inner_join(frame_gene, samp_table, by='SUBJID')

s = test_t %>%
  mutate(age = ifelse(AGE > 55, 'old', 'young')) %>%
  group_by(genotype, age) %>%
  dplyr::summarise(mu_g = mean(exp))
s$gene = gene

s

data_all <- "Y:/analysis/JSD_eQTL_new/All_Tissue_lm_results_new.tsv"

frame2 <- vroom(data_all)

gt50 = c('Whole_Blood','Stomach',
         'Colon_Sigmoid',
         'Esophagus_Gastroesophageal_Junction',
         'Colon_Transverse',
         'Artery_Aorta',
         'Heart_Atrial_Appendage',
         'Breast_Mammary_Tissue',
         'Prostate',
         'Heart_Left_Ventricle',
         'Cells_Culturedfibroblasts',
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

# data_best6 <- "Y:/analysis/JSD_eQTL/best6_lm_results.tsv"
# frame <- vroom(data_best6)

frame2$Tissue %>% unique()

tissue_names <-  vroom('../analysis/misc/selected_tissue_abbreviations_etl.tsv')

frame2$Tissue_ab <- mapvalues(frame2$Tissue, tissue_names$SMTSD, tissue_names$tissue)

g = ggplot(frame2 %>% 
             filter(age == 'neut') %>%
             group_by(Tissue_ab) %>%
             mutate(mu = mean(slope_AGE))
             , aes(x = reorder(Tissue_ab, mu), y = slope_AGE))

g0 = g + stat_summary()

ggsave(paste0(path_temp, "slope_age.png"),  plot=g0, width = 20, height = 12, scale=1, dpi=600)

t_sum = frame2 %>% 
  mutate(P=pval_combined) %>%
  group_by(Tissue,age) %>%
  arrange(P) %>%
  dplyr::mutate(exp_P=row_number()/n(),
         mu_P = mean(P)) %>%
  ungroup() %>%
  dplyr::select(P,exp_P,gene,Tissue,age) %>%
  mutate(logP=-log(P)) %>%
  arrange(age)

t_sum_gene = t_sum %>% 
  dplyr::select(age, gene, Tissue, logP) %>%
  unique() %>%
  pivot_wider(names_from = age, values_from = logP, values_fn=mean) %>%
  filter(!is.na(young)) %>%
  mutate(delta = old - young) %>%
  arrange(delta)

t_sum_gene %>%
  mutate(up = ifelse(delta < 0, 'up', 'down')) %>%
  group_by(up) %>%
  dplyr::summarise(n = n())

t_sum_blood = t_sum_gene %>% filter(Tissue == "Whole_Blood")

t_sum_blood$gene

t_summary_stats = frame2 %>% 
  dplyr::mutate(P=pval_combined) %>%
  dplyr::mutate(logP = -log(P, 10)) %>%
  dplyr::group_by(Tissue,age) %>%
  dplyr::summarize(mu_logP = mean(logP,na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Tissue,age,mu_logP) %>%
  filter(!is.nan(mu_logP)) %>%
  pivot_wider(names_from="age",values_from="mu_logP") %>%
  mutate(delta=young-old) %>%
  filter(Tissue %in% gt50) %>%
  arrange(-delta)
  
t_summary_stats_slope = frame2 %>% 
    dplyr::mutate(slope = `slope_genotype[T.0|1]`) %>%
    dplyr::group_by(Tissue,age) %>%
    dplyr::summarize(mu_slope = mean(slope, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(Tissue,age, mu_slope) %>%
    filter(!is.nan(mu_slope)) %>%
    pivot_wider(names_from="age",values_from="mu_slope") %>%
    mutate(delta=young-old) %>%
    filter(Tissue %in% gt50) %>%
    arrange(-delta)

t_summary_new = inner_join(t_summary_stats, t_summary_stats_slope, by='Tissue')

t_dist = frame2 %>% 
  filter(Tissue %in% gt50) %>%
  dplyr::mutate(slope = `slope_genotype[T.0|1]`) %>%
  dplyr::mutate(logP = -log(pval_combined, 10)) %>%
  dplyr::select(Tissue, age, slope, gene, logP) %>%
  filter(!is.nan(slope)) %>%
  pivot_wider(names_from = "age", values_from = c("slope", "logP"), values_fn = mean ) %>%
  mutate(delta_slope = slope_young - slope_old) %>%
  mutate(delta_P = logP_young - logP_old) 

g = ggplot(t_dist, aes(x = delta_P, y = delta_slope, col = Tissue))
g + stat_summary_bin() + geom_smooth(method = 'lm', se=FALSE) + ylim(c(-0.5, 0.5))

m0 = t_dist %>% group_by(Tissue) %>%
  do(tidy(lm(delta_slope~delta_P, data=.)))

m0 %>% filter(estimate < 0)

a = m0 %>% filter(p.value < 0.05 / 27) 

a$Tissue

g = ggplot(t_summary_new)

R = cor.test(t_summary_new$delta.x, t_summary_new$delta.y)

R

g0 = g + geom_point(aes(x = delta.x, y = delta.y, fill = Tissue), size=5, pch=21, colour='black') +
  geom_text_repel(aes(x = delta.x, y = delta.y, label = Tissue), force = 1, nudge_y = 0.001, size=5) + 
  labs(
    x=TeX('$\\Delta{log P value}'),
    y=TeX('$\\Delta{slope}$')
  )+
  theme_bw(base_size=16)+
  theme(
    legend.position = "None",
    legend.box = "vertical",
    text = element_text(size=14),
    axis.title = element_text(size=16),
    axis.text = element_text(size=16)
  ) + 
  geom_text(aes(x = -0.6, y = 0.0055, label = paste('R =', signif(R$estimate, 2), 'P = ', signif(R$p.value, 2) ) ), size=6)

g0

ggsave(paste0(path_temp, "JSD_logP.png"),  plot=g3, width = 14, height = 12, scale=0.5, dpi=600)



t_summary_stats$Tissue_ab <- mapvalues(t_summary_stats$Tissue, tissue_names$SMTSD, tissue_names$tissue)

t_plot = inner_join(t_sum,t_summary_stats,by="Tissue")

t_plot$Tissue = factor(t_plot$Tissue, levels=t_summary_stats$Tissue)

t_plot$Tissue_ab <- mapvalues(t_plot$Tissue, tissue_names$SMTSD, tissue_names$tissue)

t_plot$Tissue_ab = factor(t_plot$Tissue_ab, levels=t_summary_stats$Tissue_ab)

t_plot$Tissue_ab


t_plot = t_plot %>%
  mutate(
    Tissue.type = Tissue %>% 
      str_split("_") %>% 
      sapply("[[", 1) %>% 
      clean.subtypes,
    Tissue.subtype = Tissue %>% 
      str_split(., pattern = "_", n=2) %>% 
      sapply(function(x){ifelse(length(x)>1, x[[2]], '')}) %>% 
      clean.subtypes,
    Tissue.neat = str_c(Tissue.type, "\n", Tissue.subtype) %>% 
      str_trim()
  )

pval_star <- function(pval) {
  if(pval > 0.001) {
    return('ns')
  } else {
    return(paste("P = ", signif(pval, 2) ) )
  }
}

a = t_plot %>% 
  dplyr::select(gene, Tissue_ab) %>%
  unique() %>%
  dplyr::group_by(Tissue_ab) %>%
  mutate(count = n()) %>%
  dplyr::select(count, Tissue_ab) %>%
  unique

a

t_plot

a

a %>% arrange(-count)

b = t_plot %>% 
  dplyr::group_by(Tissue_ab,delta) %>%
  mutate(logP = -log(P,10)) %>%
  do(w=wilcox.test(logP~age,data=.,paired=FALSE))



c = b$w[[1]]

c$statistic[1]

t_label = t_plot %>% 
  dplyr::group_by(Tissue_ab,delta) %>%
  mutate(logP = -log(P,10)) %>%
  filter(age %in% c('young', 'old')) %>%
  do(w=wilcox.test(logP~age,data=.,paired=FALSE)) %>%
  dplyr::summarize(Tissue_ab,delta, W=w$p.value, stat=w$statistic[1] )%>%
  rowwise() %>%
  mutate(p_label = pval_star(W )) %>%
  mutate(log_W = -log(W, 10)) %>%
  mutate(up_down = ifelse(delta > 0, 1, -1)) %>%
  mutate(dir_W = up_down*log_W) %>%
  arrange(-dir_W)

t_plot$Tissue_ab = factor(t_plot$Tissue_ab, levels=t_label$Tissue_ab)

t_plot$Tissue_ab

t_plot$Tissue %>% unique()

g=ggplot(t_plot %>%
           filter(!Tissue %in% c('Whole_Blood', 'Stomach')) %>%
           filter(age %in% c('young', 'old'))
         )

g2 = g+geom_point(aes(x=-log(exp_P,10),y=logP,color=age),size=.1)+
  theme_bw(base_size = 14)+
  facet_wrap(~Tissue_ab)+
  geom_text(aes(x=1.5,y=500,label=p_label),
            data=t_label %>%
              filter(!Tissue_ab %in% c('Blood', 'Stomach'))
            ,size=4) + 
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text = element_text(size=10))+
  labs(x = '-log(Expected)', y='-log(P values)') + 
  scale_color_brewer(palette="Set1")

g2

ggsave(paste0(path_temp, "quantile_eQTL.png"),  plot=g2, width = 14, height = 12, scale=0.7, dpi=600)

g = ggplot(t_plot %>% 
             filter(Tissue %in% c('Whole_Blood', 'Stomach')) %>%
             filter(age %in% c('young', 'old'))
           
           )

g2.2 = g+geom_point(aes(x=-log(exp_P,10),y=logP,color=age),size=1)+
  theme_bw(base_size = 14)+
  facet_wrap(~Tissue_ab)+
  geom_text(aes(x=1.5,y=500,label=p_label),data=t_label %>% filter(Tissue_ab %in% c('Blood', 'Stomach')),size=8) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size=16))+
  labs(x = '-log(Expected)', y='-log(P values)') + 
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_brewer(palette="Set1")

g2.2

ggsave(paste0(path_temp, "quantile_eQTL_WB.png"),  plot=g2.2, width = 14, height = 12, scale=0.5, dpi=600)

load("joint_frame.Rda")

frame_g2 = inner_join(joint.frame, t_summary_stats, by='Tissue_ab')

joint.frame$Tissue.neat

b = inner_join(joint.frame, a, by='Tissue.neat')

c = b %>% dplyr::select(pal.tissue, Tissue)

write.csv(c, '../analysis/pal_tissue.csv')

save(frame_g2, file="frame_g2.Rda")

g = ggplot(frame_g2)

pal.tissue <- readRDS('../analysis/pal_tissue.Rdata')

a = data.frame(pal.tissue)

library(data.table)
library(vroom)

fwrite(a, '../analysis/pal_tissue.csv', eol="\r\n")

vroom_write

a$Tissue.neat = rownames(a)

write.table(a, '../analysis/pal_tissue.csv')

txt = cor.test(frame_g2$delta, frame_g2$diff)

txt

g3 = g + geom_point(aes(x = delta, y = diff, fill = Tissue.neat), size=5, pch=21, colour='black') +
  geom_text_repel(aes(x = delta, y = diff, label = Tissue_ab), force = 1, nudge_y = 0.001, size=5) + 
  labs(
    x=TeX('$\\Delta{log P value}'),
    y=TeX('$\\Delta{JSD}$')
  )+
  theme_bw(base_size=16)+
  theme(
    legend.position = "None",
    legend.box = "vertical",
    text = element_text(size=14),
    axis.title = element_text(size=16),
    axis.text = element_text(size=16)
  ) + 
  scale_fill_manual(
    values = pal.tissue
  ) +
  geom_text(aes(x = -0.6, y = 0.014, label = 'R = 0.523, P = 0.005'), size=6)

g3

ggsave(paste0(path_temp, "JSD_logP.png"),  plot=g3, width = 14, height = 12, scale=0.5, dpi=600)


# 
# gene_lst = frame2 %>% 
#   filter(Tissue %in% c('Whole_Blood')) %>%
#   dplyr::select(R2,Tissue,age,gene) %>%
#   pivot_wider(names_from="age", values_from ="R2") %>%
#   mutate(delta=young-old) %>%
#   arrange(-delta) 
# 
# organism = "org.Hs.eg.db"
# library(organism, character.only = TRUE)
# organism = org.Hs.eg.db
# 
# gene_lst = gene_lst %>%
#   rowwise() %>%
#   mutate(gene = strsplit(gene, "\\."))
# gene_lst = gene_lst %>%
#   rowwise() %>% 
#   mutate(gene = gene[1])
# 
# genes_wb = gene_lst$delta
# 
# names(genes_wb) = gene_lst$gene
#   
# gs_wb <-  gseGO(geneList=genes_wb, 
#                   ont ="ALL", 
#                   keyType = "ENSEMBL", 
#                   nPerm = 10000, 
#                   minGSSize = 3, 
#                   maxGSSize = 800, 
#                   pvalueCutoff = 0.05, 
#                   verbose = TRUE, 
#                   OrgDb = organism, 
#                   pAdjustMethod = "none")
#   
# dotplot(gs_wb, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("Whole Blood")
# 
# library(PNWColors)
# 
# pal=pnw_palette("Bay",100)
# 
# d = data.frame(gs_wb)
# 
# head(d$core_enrichment)
# 
# d$count = strsplit(d$core_enrichment, '/')
# 
# d$count_size = lapply(d$count, length)
# 
# 
# d$count = unlist(d$count_size)
# 
# d = d %>% mutate(perc = count / setSize)
# 
# d_top = d %>% filter(p.adjust < 0.05 & enrichmentScore > 0.0) %>% top_n(10, -p.adjust)
# 
# d_bot = d %>% filter(p.adjust < 0.05 & enrichmentScore < 0.0) %>% top_n(10, -p.adjust)
# 
# g=ggplot(rbind(d_top, d_bot) %>% arrange(NES))
# 
# g2 = g + geom_segment(aes(x=0,
#                      xend=-enrichmentScore,
#                      y=reorder(Description,enrichmentScore),
#                      yend=reorder(Description,enrichmentScore)),
#                  size=.3,
#                  color="grey") +
#   geom_point(aes(x=0,
#                  y=reorder(Description,enrichmentScore)),
#              size=.3,
#              color="grey") +
#   geom_point(aes(x=-enrichmentScore,
#                  y=reorder(Description,enrichmentScore),
#                  fill=pvalue,
#                  size=setSize),
#              shape=24,
#              stroke=.2)+
#   theme_bw(base_size=16)+
#   scale_fill_gradientn(colors=pal)+
#   scale_y_discrete("")+
#   scale_size(range=c(2,8))+
#   guides(size=FALSE)+
#   theme(legend.position=c(0.82,0.85),
#         legend.key.size=unit(.3,"cm"),
#         legend.box.background = element_rect(colour = "grey")) +
#   ggtitle('Whole Blood') + 
#   labs(x = 'Enrichment score')
# 
# ggsave(paste0(path_temp, "WB_gsea_eQTL.png"),  plot=g2, width = 16, height = 12)
# 


















# 
# g1 = ggplot(frame %>% mutate(log_p = -log10(`p_val_genotype[T.1|0]`)) ) +
#   geom_point(aes(x = `slope_genotype[T.0|1]`, y = log_p, col=age)) +
#   # geom_text_repel(
#   #   data = subset(frame2 %>% mutate(log_p = -log10(`p_val_genotype[T.1|0]`)), pval_combined < 1e-250),
#   #   aes(x = slope_combined, y = log_p, label = gene),
#   #   size = 4,
#   #   box.padding = unit(0.35, "lines"),
#   #   point.padding = unit(0.3, "lines")
#   # ) +
#   scale_color_manual(values=pnw_palette("Sunset2",2)) + 
#   facet_wrap(~Tissue) +
#   theme_bw()
# 
# g1
# 
# t = frame2
# 
# t_sum1 = t %>% mutate(P=pval_combined) %>%
#   arrange(P) %>%
#   group_by(age, Tissue) %>%
#   mutate(exp_P=row_number()/n()) %>%
#   dplyr::select(P,exp_P,gene) %>%
#   mutate(logP=-log10(P)) 
# 
# g=ggplot(t_sum1)
# 
# g3 = g+geom_point(aes(y=-log(P),x=-log(exp_P), color=age),alpha=.5)+
#   theme_bw()+
#   facet_wrap(~Tissue) + 
#   geom_text() + 
#   scale_color_brewer(palette = "Set1") 
# 
# g3
# 
# frame2 = frame %>%
#   rowwise() %>%
#   mutate(pval_combined = fisher(`p_val_genotype[T.0|1]`, `p_val_genotype[T.1|0]`, `p_val_genotype[T.1|1]`) ) %>%
#   mutate(slope_combined = (`slope_genotype[T.0|1]` + `slope_genotype[T.1|0]` + `slope_genotype[T.1|1]`) / 3 )
# 
# write.csv(frame2, "Y:/analysis/JSD_eQTL/best6_lm_results_2.tsv")
# 
# frame2$pval_combined
# 
# frame2$Significant <- ifelse(frame2$pval_combined < 1e-5, "Sig", "Not Sig")
# 
# path_temp <- 'G:/Shared drives/sudmantlab/projects/agingGeneRegulation/figures/figure1/fig1_'
# 
# g1 = ggplot(frame2 %>% mutate(log_p = -log10(pval_combined)) ) +
#   geom_point(aes(x = slope_combined, y = log_p, col=age)) +
#   geom_text_repel(
#     data = subset(frame2 %>% mutate(log_p = -log10(pval_combined)), pval_combined < 1e-250),
#     aes(x = slope_combined, y = log_p, label = gene),
#     size = 4,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines")
#   ) +
#   scale_color_manual(values=pnw_palette("Sunset2",2)) + 
#   facet_wrap(~Tissue) +
#   theme_bw()
# 
# g1
# 
# ggsave(paste0(path_temp, "volcano_eQTL.png"),  plot=g1, width = 14, height = 12, scale=0.5)
# 
# t = frame2
# 
# t_sum1 = t %>% mutate(P=pval_combined) %>%
#   arrange(P) %>%
#   group_by(age, Tissue) %>%
#   mutate(exp_P=row_number()/n()) %>%
#   dplyr::select(P,exp_P,gene) %>%
#   mutate(logP=-log10(P)) 
# 
# g=ggplot(t_sum1)
# 
# g2 = g+geom_point(aes(y=-log(P),x=-log(exp_P), color=interaction(age, Tissue) ),alpha=.5)+
#   theme_bw()+
#   scale_color_brewer(palette = "Set1") 
# 
# g2
# 
# ggsave(paste0(path_temp, "QQ_eQTL.png"),  plot=g2, width = 14, height = 12, scale=0.5)
# 
# g3 = g+geom_point(aes(y=-log(P),x=-log(exp_P), color=age),alpha=.5)+
#   theme_bw()+
#   facet_wrap(~Tissue) + 
#   geom_text() + 
#   scale_color_brewer(palette = "Set1") 
# 
# g3
# 
# ggsave(paste0(path_temp, "QQ_eQTL.png"),  plot=g3, width = 14, height = 12, scale=0.5)
# 
# t_test = frame2 %>% group_by(Tissue) %>%
#   do(tidy(t.test(pval_combined ~ age, data=.)))
# 
# 
# 
# ggplot(frame %>% mutate(log_p = -log10(pval_combined))  %>% filter(age == 'old') ) +
#   geom_point(aes(x = slope_combined, y = log_p, col=Significant)) +
#   geom_text_repel(
#     data = subset(frame %>% mutate(log_p = -log10(pval_combined))  %>% filter(age == 'old'), pval_combined < 0.00001),
#     aes(x = slope_combined, y = log_p, label = gene),
#     size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines")
#   ) + theme_bw() + ggtitle('old volcano')
# 
# ggplot()
# 
# ggplot(frame %>% mutate(log_p = -log10(pval_combined)))  +
#   geom_boxplot(aes(y = log_p, x=age)) +
#   stat_compare_means(aes(y = log_p, x=age), method = "t.test")
# test_t = frame %>%
# mutate(genotype = chr17_40050862_C_G_b38, exp = ENSG00000108342.12, SUBJID = ...1)
# 
# test_t = inner_join(test_t, samp_table, by='SUBJID')
# 
# test_t %>%
#   group_by(genotype) %>%
#   dplyr::summarize(count = n())
# 
# g1 = ggplot(test_t %>%
#               filter(!is.na(exp)) %>%
#               mutate(age = ifelse(AGE > 55, 'old', 'young')) %>%
#               mutate(genotype = ifelse(genotype == '0|1', '1/0', genotype)) %>%
#               mutate(genotype = ifelse(genotype == '1|0', '1/0', genotype)) %>%
#               mutate(genotype = ifelse(genotype == '0|0', '0/0', genotype)) %>%
#               mutate(genotype = ifelse(genotype == '1|1', '1/1', genotype))
# ) +
#   stat_summary(aes(y = exp, x = fct_rev(age), col= age), size = 2, alpha = 0.8) +
#   geom_jitter(aes(y = exp, x = fct_rev(age), col= age) , alpha = 0.3) +
#   facet_wrap(~genotype, nrow = 1) + 
#   labs(
#     x=TeX('genotype'),
#     y=TeX('expression')
#   )+
#   theme_bw(base_size=16)+
#   theme(
#     legend.position = "None",
#     legend.box = "vertical",
#     text = element_text(size=14),
#     axis.title = element_text(size=16),
#     axis.text = element_text(size=16),
#     strip.background = element_blank(),
#     strip.text = element_text(size=16)
#   ) + 
#   scale_color_brewer(
#     palette = 'Set1'
#   )
# 
# g1
# 
# ggsave(paste0(path_temp, 'gene_example.png'), plot = g1, width = 14, height = 12, scale=0.5, dpi=600)
# ggplot(frame %>% mutate(log_p = -log10(pval_combined)) )  + stat_ecdf(aes(x = log_p, col=age))
