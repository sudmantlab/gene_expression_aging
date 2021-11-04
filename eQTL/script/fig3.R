library(datasets)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)
library(vroom)
library(tidyverse)
library(broom)
library(latex2exp)
library(ggtern)
library(viridis)
library(gplots)
source('../Rscript/cleanSubType.R')

path = "Y:/analysis/Predixcan_eQTL/"

path_temp <- 'G:/Shared drives/sudmantlab/projects/agingGeneRegulation/figures/figure3/fig3_'

frame = vroom(paste0(path, "r2_b27.tsv"))

pval_star <- function(pval) {
  if(pval > 0.001) {
    return('ns')
  } else {
    return(paste("P = ", signif(pval, 2) ) )
  }
}

b = frame  %>% 
  dplyr::group_by(Tissue) %>%
  do(tidy(t.test(R2~age,data=., paired=TRUE)) ) %>%
  arrange(estimate)

frame$Tissue = factor(frame$Tissue, levels=b$Tissue)

frame$Tissue

# g1 = ggplot(frame)
# 
# g1 + geom_violin(aes(x=age, y=R2, fill=age)) +
#   facet_wrap(~Tissue) +
#   theme_classic() +
#   ylim(c(-0.5, 1)) + 
#   geom_text(aes(x = 1.5, y = 0.7, label = p.value), size=3, data = b %>% mutate(p.value = pval_star(p.value))) + 
#   geom_text(aes(x = 1.5, y = 0.87, label = paste("delta = ", signif(-1*estimate, 2) )), size = 3, data= b) + 
#   theme(strip.text = element_text(size=7))
# 
# ggsave(paste0(path_temp, "facet_dist_b27.pdf"))


load("Z:/sudmant_lab/gene_expression_variance/Rscript/frame_g2.Rda")

a = inner_join(frame_g2 %>% 
                 mutate(Tissue = Tissue.y) %>%
                 dplyr::select(Tissue, delta, young, old, Tissue_ab, Tissue.neat, estimate),
               b %>% mutate(estimate_n = estimate), by = "Tissue")

a 

cor.test(a$delta, a$estimate_n)

cor.test(a$estimate.x, a$estimate_n)

pal.tissue <- readRDS('../analysis/pal_tissue.Rdata')

ggplot(a) + 
  geom_point(aes(x = delta, y = estimate_n, fill = Tissue.neat), size=5, pch=21, colour='black') +
  geom_text_repel(aes(x = delta, y = estimate_n, label = Tissue_ab), force = 1, nudge_y = 0.001, size=5) + 
  labs(
    x=TeX('$\\Delta{log P value}'),
    y=TeX('$\\Delta{h^2}$')
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
  geom_text(aes(x = -0.6, y = 0.006, label = 'R = 0.569, P = 0.00193'), size=6)

ggsave(paste0(path_temp, "Predixcan_vs_eQTL_b27.pdf"))
tissue_names <-  vroom('../analysis/misc/selected_tissue_abbreviations_etl.tsv')
a$Tissue_ab <- mapvalues(a$Tissue, tissue_names$SMTSD, tissue_names$tissue)

ggplot(a %>% mutate(
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
)) + 
  geom_point(aes(x = estimate.x, y = estimate.y, fill = Tissue.neat), size=5, pch=21, colour='black') +
  geom_text_repel(aes(x = estimate.x, y = estimate.y, label = Tissue_ab), force = 1, nudge_y = 0.001, size=5) + 
  labs(
    x=TeX('$\\Delta{JSD} (old - young)'),
    y=TeX('$\\Delta{h^2}$ (old - young )')
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
  geom_text(aes(x = 0, y = 0.006, label = 'R = -0.597, P = 9.9e-4'), size=6)

#ggsave(paste0(path_temp, "Predixcan_vs_JSD_b27.pdf"))

cor.test(a$estimate.x, a$delta)

ggplot(a) + 
  geom_point(aes(x = estimate.x, y = -delta, fill = Tissue.neat), size=5, pch=21, colour='black') +
  geom_text_repel(aes(x = estimate.x, y = -delta, label = Tissue_ab), force = 1, nudge_y = 0.001, size=5) + 
  labs(
    x=TeX('$\\Delta{JSD}'),
    y=TeX('$\\Delta{log P value}$')
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
  geom_text(aes(x = 0.0004, y = 0.65, label = 'R = -0.48, P = 0.01'), size=6) +
  coord_flip()

#ggsave(paste0(path_temp, "eQTL_vs_JSD_b27.pdf"), scale = 0.8, width = 10, height = 8)

ggplot(a %>% 
         mutate(p_star = ifelse(p.value < 0.001, '**', 'ns'))
       , 
       aes(x = estimate_n,
           y = reorder(Tissue, estimate_n),
           fill = Tissue.neat)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(
    values = pal.tissue
  ) +
  theme_classic() +
  labs(y = 'Tissue', x = 'delta R2 (young - old)') + 
  geom_text(aes(x = -0.0055, y = Tissue, label = p_star), size = 6) + 
  theme(legend.position = 'None')

#ggsave(paste0(path_temp, "summary_b27.pdf"))

t = vroom("G:/Shared drives/sudmantlab/projects/agingGeneRegulation/tables/agecorr_res.csv")

names(t)

t_dat = t %>% 
  mutate(G =`Genetic R^2`, 
         A = `Age R^2`,
         E = 1-G-A,
         Coef = `Age Coef`
         ) %>%
  mutate(G= ifelse(G<0,0,G)) %>%
  mutate(A= ifelse(A<0,0,A)) %>%
  dplyr::select(G,A,E,Gene,Tissue, gene, Coef) 

t_sum = t_dat %>%
  filter(E<1) %>%
  pivot_longer(cols = c(G,A,E), names_to="var_component")

t_sum = t_dat %>% 
  #filter(E<1) %>%
  group_by(Tissue) %>%
  mutate(A_mu = mean(A),
         G_mu = mean(G),
         E_mu = mean(E)) %>%
  ungroup() %>%
  pivot_longer(cols = c(G,A,E), names_to="var_component") %>%
  #filter(var_component!="E") %>%
  group_by(var_component,Tissue,A_mu,G_mu,E_mu) %>%
  dplyr::summarize(mu = mean(value),
            med = median(value),
            sigma = sd(value),
            count = n()) %>%
  mutate(se = sigma/sqrt(count))

a = inner_join(frame_g2 %>% mutate(Tissue = Tissue.y),
               t_sum %>% filter(var_component == "A") %>% mutate(Tissue = Tissue),
               by = 'Tissue')

cor.test(a$A_mu, a$estimate)

ggplot(a)+ 
  geom_point(aes(x = estimate, y = A_mu, fill = Tissue.neat), size=5, pch=21, colour='black') +
  geom_text_repel(aes(x = estimate, y = A_mu, label = Tissue_ab), force = 1, nudge_y = 0.001, size=5) + 
  labs(
    x=TeX('$\\Delta{JSD}'),
    y=TeX('$\\R^2(Age)$')
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
  geom_text(aes(x = -0.00035, y = 0.07, label = 'R = 0.457, P = 0.0165'), size=6)

#ggsave(paste0(path_temp, "AGER2_vs_JSD_b27.pdf"))

# g=ggplot(t_sum)
# g+geom_point(aes(y=mu,x=reorder(tissue,A_mu),color=var_component),size=2)+
#   geom_errorbar(aes(ymin=mu-se,ymax=mu+se,x=reorder(tissue,mu),color=var_component),width=0)+
#   #facet_wrap(~var_component,ncol=1)+
#   theme_bw()+
#   coord_flip()+
#   scale_color_brewer(palette="Set1") + 
#   labs(x = "tissues", y = "mean variance explained", color = "variable")

t_sum   %>% arrange(Tissue)
g=ggplot(t_sum %>% filter(var_component!="E") %>% mutate(Tissue = Tissue,
  Tissue.type = Tissue %>% 
    str_split("_") %>% 
    sapply("[[", 1) %>% 
    clean.subtypes,
  Tissue.subtype = Tissue %>% 
    str_split(., pattern = "_", n=2) %>% 
    sapply(function(x){ifelse(length(x)>1, x[[2]], '')}) %>% 
    clean.subtypes,
  Tissue.neat = str_c(Tissue.type, " ", Tissue.subtype) %>% 
    str_trim()
) )
g+geom_point(aes(y=mu,x=reorder(Tissue.neat,A_mu),color=var_component), size=4)+
  geom_errorbar(aes(ymin=mu-se,ymax=mu+se,x=reorder(Tissue.neat,mu),color=var_component), width=0.7, color = 'black')+
  #facet_wrap(~var_component,ncol=1)+
  theme_bw()+
  coord_flip()+
  scale_color_brewer(palette="Set2", labels = c('Age', 'Genetics'))+ 
  labs(x = "tissues", y = TeX('$\\mu_{R^2}$'), color = TeX('$R^2$')) + 
  theme(axis.title=element_text(size=14,face="bold"))

ggsave(paste0(path_temp, "summary_A_G_b27.pdf"), scale = 0.7)

t_tern_sum  = t_dat  %>% 
  group_by(Gene) %>%
  summarise(A=mean(A),
            G=mean(G),
            `E+O`=mean(E))

g1 = ggtern(data=t_dat %>% filter(E < 1) %>%
         mutate(`E+O` = E) %>%
         mutate(`A+G` = A+G)
         )+
  geom_point(aes(A,G,`E+O`,col = `A+G`),alpha=0.3,shape=21,size=4)+
  geom_mask()+
  theme_bw(base_size=12)+
  theme_showarrows()+
  #theme(legend.position="bottom")+
  labs(Tarrow="Genetics",
       Rarrow="Environment + Others",
       Larrow="Age",
       T = TeX('$R^2_{genetics}$'),
       R = TeX('$R^2_{environment}$'),
       L = TeX('$R^2_{age}$'),
       col = TeX('$R^2_{genetics} + R^2_{age}$'),
       )+
  scale_color_viridis()

g1

ggsave(paste0(path_temp, "triangle_A_G_E_b27.png"), scale = 0.85, plot = g1, dpi=600)

summary(lm(G~A,data=t_dat %>% filter(E<1)))

t_genes = t_dat %>% filter(E < 1) %>% pivot_wider(names_from = Tissue, values_from = c(A, G, Coef))

colnames(t_genes)

t_genes[is.na(t_genes)] <- 0

t_genes_A = t_dat %>% 
  filter(E < 1) %>%
  dplyr::select(A, Tissue, Gene, gene) %>%
   pivot_wider(names_from = Tissue, values_from = c(A))

t_genes_Coef = t_dat %>% 
  filter(E < 1) %>%
  dplyr::select(Coef, Tissue, Gene, gene) %>%
  pivot_wider(names_from = Tissue, values_from = c(Coef))

#t_genes_A[is.na(t_genes_A)] <- 0

t_genes_G = t_dat %>% filter(E < 1) %>% 
  dplyr::select(G, Tissue, Gene, gene) %>%
  pivot_wider(names_from = Tissue, values_from = c(G))

t_genes_G[is.na(t_genes_G)] <- 0

b = t_genes_A %>% drop_na() %>% 
  mutate(mu = rowMeans(dplyr::select(., -c(Gene, gene)) ))

t_genes_A

b = b %>%
  rowwise() %>%
  mutate(gene = str_remove(Gene, "\\..*")) %>%
  filter(! gene == "")

c = c %>%
  rowwise() %>%
  mutate(gene = str_remove(Gene, "\\..*")) %>%
  filter(! gene == "")


d = t_genes_Coef %>% drop_na() %>%
  mutate(mu = rowMeans(dplyr::select(., -c(Gene, gene)) )) %>%
  rowwise() %>%
  mutate(gene = str_remove(Gene, "\\..*")) 
# 
# long_gene_A = b %>%
#   pivot_longer(!c(Gene, gene, gene, mu), "names_to" = "Tissue") %>%
#   group_by(Tissue) %>%
#   mutate(mu_T = mean(value))
# 
# g2 = ggplot(long_gene_A %>% top_frac(0.5, wt=mu)) + 
#   geom_tile(aes(x = reorder(Tissue, mu_T), y = reorder(gene, mu), fill = value)) +
#   scale_x_discrete(labels = abbreviate) +
#   scale_fill_viridis()
# 
# g2
# 
# ggsave(paste0(path_temp, "heatmap_A_b27.png"), plot = g2, dpi=600, height = 24, width = 14)

library(matrixStats)

tis_sp_score = function(x) {
  ma = max(x)
  sc = 0
  for (i in x) {
    sc = sc + (1 - i / ma)
  }
  sc = sc / length(x)
  return(sc)
}

colfunc <- colorRampPalette(c("white", "red"))

mat_col = data.matrix(b %>% dplyr::select(-c(Gene, gene, mu)))

tsp_A = data.frame(apply(mat_col, 1, tis_sp_score) )

vars_A = data.frame(rowVars(mat_col))

tsp_A$gene = b$gene

vars_A$gene = b$gene

# rownames(mat_col) = b$gene
# 
# colnames(mat_col) = abbreviate(colnames(mat_col))
# 
# b$gene


# png(paste0(path_temp, "heatmap_A_b27_cluster_new.png"), height = 7*300, width = 5*300, res = 300, pointsize = 8)
# heatmap.2(mat_col, scale = "none", srtCol = 45, col = colfunc(100), dendrogram = "row",
#           hclustfun = function(x) hclust(x, method="average"),
#           keysize = 0.7, key.title = '',
#           key.xlab = 'Age associated R2',
#           trace = "none", density.info = "none")
# dev.off()


c = t_genes_G %>% drop_na() %>% 
  mutate(mu = rowMeans(dplyr::select(., -c(Gene, gene)) )) %>%
  filter(! gene == "") %>%
  rowwise() %>%
  mutate(gene = str_remove(Gene, "\\..*")) %>%
  filter(! gene == "")

mat_col_G = data.matrix(c %>% dplyr::select(-c(Gene, gene, mu)))

tsp_G = data.frame(apply(mat_col_G, 1, tis_sp_score) )

vars_G = data.frame(rowVars(mat_col_G))

vars_G$gene = c$gene

tsp_G$gene = c$gene

rownames(mat_col_G) = c$gene

colnames(mat_col_G) = abbreviate(colnames(mat_col_G))

# png(paste0(path_temp, "heatmap_G_b27_cluster_new.png"), height = 7*300, width = 5*300, res = 300, pointsize = 8)
# heatmap.2(mat_col_G, scale = "none", col = colfunc(100), srtCol = 45, dendrogram = "row",
#           hclustfun = function(x) hclust(x, method="average"),
#           keysize = 0.7, key.title = '',
#           key.xlab = 'Genetic associated R2',
#           trace = "none", density.info = "none")
# dev.off()

vars_AG = inner_join(vars_A, vars_G, by = 'gene') 

vars_AG = vars_AG %>%
  mutate(A = rowVars.mat_col., G = rowVars.mat_col_G.) %>%
  select(gene, A, G)

t.test(vars_AG$A, vars_AG$G, paired = TRUE, alternative = 'two.sided')

tsp_AG = inner_join(tsp_A, tsp_G, by = 'gene')

tsp_AG = tsp_AG %>%
  mutate(A = apply.mat_col..1..tis_sp_score., G = apply.mat_col_G..1..tis_sp_score.) %>%
  dplyr::select(gene, A, G)

t.test(tsp_AG$A, tsp_AG$G, paired = TRUE, alternative = 'two.sided')

library(ggpubr)

g3 = ggplot(tsp_AG %>%
         pivot_longer(!gene, names_to = 'name', values_to='value') %>%
         filter(value != 0)) +
  geom_boxplot(aes(x = name, y = value, fill = name)) +
  # stat_summary(aes(x = name, y = value, col = name), alpha = 0.8, size = 5) + 
  # geom_jitter(aes(x = name, y = value, col = name), alpha = 0.3) + 
  scale_fill_brewer(palette = 'Set2') + 
  #coord_cartesian(ylim = c(0, 0.01)) + 
  theme_classic(base_size = 14) + 
  theme(legend.position = 'None') + 
  labs(x = TeX('$R^2$'), y = "tissue specificity score", fill = "attribute") +
  geom_text(aes(x = 1.5, y = 0.03), label = "p < 2.2e-16", size = 5) +
  scale_x_discrete(labels = c('Age', 'Genetics'))

g3

ggsave(paste0(path_temp, "_row_tsp.png"), plot = g3, width = 14, height = 12, dpi = 300, scale = 0.5)

# 

# b = t_genes_A[which(rowMeans(!is.na(t_genes_A)) > 0.9), ]%>%
#   rowwise() %>%
#   mutate(gene = str_remove(Gene, "\\..*"))


library(PNWColors)

pal=pnw_palette("Bay",100)

# g = ggplot(gs_res %>% top_n(30, wt=NES))
# 
# g +  geom_point(aes(x=NES,
#                  y=reorder(Description,NES),
#                  fill=p.adjust,
#                  size=setSize),
#              shape=24,
#              stroke=.2)+
#   theme_bw(base_size=16)+
#   scale_fill_gradientn(colors=pal)+
#   scale_y_discrete("")+
#   scale_size(range=c(2,8))+
#   theme(legend.key.size=unit(.3,"cm"),
#         legend.box.background = element_rect(colour = "grey")) +
#   labs(x = 'Normalized Enrichment score')
# 
# ggsave(paste0(path_temp, "GO_pLI_lt_0.25_T5.pdf"))
# 
# c = t_genes_A[which(rowMeans(!is.na(t_genes_A)) > 0.9), ]%>%
#   replace(is.na(.), 0) %>%
#   dplyr::mutate(mu = rowMeans(dplyr::select(., -c(Gene, gene)) ))
# 
# c = b %>%
#   dplr::mutate(mu = rowMeans(dplyr::select(., -c(Gene, gene)) ))

b  = b %>%
  rowwise() %>%
  dplyr::mutate(gene = str_remove(Gene, "\\..*"))
#
# b = t_genes_A[which(rowMeans(!is.na(t_genes_A)) > 0.9), ]%>%
#   rowwise() %>%
#   mutate(gene = str_remove(Gene, "\\..*"))

genelist = d$mu

names(genelist) = d$gene

genelist = sort(genelist, decreasing = TRUE)


library(clusterProfiler)
library(enrichplot)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
organism = org.Hs.eg.db


gs <-  gseGO(geneList=genelist,
             ont ="BP",
             keyType = "ENSEMBL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.025,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "fdr")

gs_res = data.frame(gs)

#write.csv2(gs_res, file = "G:/Shared drives/sudmantlab/projects/agingGeneRegulation/tables/gsea_result_beta_b27.csv")

g = ggplot(gs_res)

g +  geom_point(aes(x=NES,
                    y=reorder(Description,NES),
                    fill=qvalues,
                    size=setSize),
                shape=24,
                stroke=.2)+
  theme_bw(base_size=16)+
  scale_fill_gradientn(colors=pal)+
  scale_y_discrete("")+
  scale_size(range=c(2,8))+
  theme(legend.key.size=unit(.3,"cm"),
        legend.box.background = element_rect(colour = "grey")) +
  labs(x = 'Normalized Enrichment score', fill = "p adjusted")

ggsave(paste0(path_temp, "GO_AGE_BP_COMMON_b27_new_0.025.pdf"), width = 16, height = 12, dpi = 300, scale = 0.8)

library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_list=names(genelist)
test=getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
           filters = "ensembl_gene_id" , values = gene_list, mart = mart)

listAttributes(mart)
test

names(gene_list) = test$entrezgene_id

genelist

names(genelist)

c = b$Gene
# 
# 
# Tissues_c7 = c('Esophagus_Mucosa', 'Colon_Transverse', 'Whole_Blood', 'Heart_Left_Ventricle',
#                'Esophagus_Gastroesophageal_Junction', 'Heart_Atrial_Appendage', 'Skin_Sun_Exposed_Lower_leg')
# 
# Tissues_c7 = c( 'Whole_Blood')
# pal=pnw_palette("Bay",100)
# 
# for (Tissue_t in Tissues_c7) {
#   
#   b = t_dat %>% 
#     filter(Tissue == Tissue_t) %>%
#     filter(!A == 0) %>%
#     dplyr::mutate(coef_sign = ifelse(Coef > 0, 1, -1)) %>%
#     dplyr::mutate(A_sign = A*coef_sign)
#   
#   genelist = b$A
#   
#   names(genelist) = b$gene
#   
#   genelist = sort(genelist, decreasing = TRUE)
#   
#   gs <-  gseGO(geneList=genelist, 
#                ont ="ALL", 
#                keyType = "ENSEMBL", 
#                minGSSize = 3, 
#                maxGSSize = 800, 
#                pvalueCutoff = 0.05, 
#                verbose = TRUE, 
#                OrgDb = organism, 
#                pAdjustMethod = "fdr")
#   
#   gs_res = data.frame(gs)
#   
#   pal=pnw_palette("Bay",100)
#   
#   g = ggplot(gs_res %>% top_n(20, wt = p.adjust))
#   
#   g_temp = g +  geom_point(aes(x=NES,
#                       y=reorder(Description,NES),
#                       fill=p.adjust,
#                       size=setSize),
#                   shape=24,
#                   stroke=.2)+
#     theme_bw(base_size=16)+
#     scale_fill_gradientn(colors=pal)+
#     scale_y_discrete("")+
#     scale_size(range=c(2,8))+
#     theme(legend.key.size=unit(.3,"cm"),
#           legend.box.background = element_rect(colour = "grey")) +
#     labs(x = 'Normalized Enrichment score')
#   
#   ggsave(paste0(path_temp, Tissue_t,"_GO_AGE_COMMON_b27.pdf"), plot = g_temp, width = 12, height = 16, dpi = 300)
#   
# }
# genelist
# a = gs_res %>% top_n(20, wt=p.adjust)
# 
# for (Tissue_t in Tissues_c7) {
#   
#   b = t_dat %>% filter(Tissue == Tissue_t) %>%
#     filter(!G == 0)
#   
#   genelist = b$G
#   
#   names(genelist) = b$gene
#   
#   genelist = sort(genelist, decreasing = TRUE)
#   
#   gs <-  gseGO(geneList=genelist, 
#                ont ="ALL", 
#                keyType = "ENSEMBL", 
#                nPerm = 10000, 
#                minGSSize = 3, 
#                maxGSSize = 800, 
#                pvalueCutoff = 1, 
#                verbose = TRUE, 
#                OrgDb = organism, 
#                pAdjustMethod = "fdr")
#   
#   gs_res = data.frame(gs)
#   
#   pal=pnw_palette("Bay",100)
#   
#   g = ggplot(gs_res %>% top_n(-20, wt = p.adjust))
#   
#   g_temp = g +  geom_point(aes(x=NES,
#                       y=reorder(Description,NES),
#                       fill=p.adjust,
#                       size=setSize),
#                   shape=24,
#                   stroke=.2)+
#     theme_bw(base_size=16)+
#     scale_fill_gradientn(colors=pal)+
#     scale_y_discrete("")+
#     scale_size(range=c(2,8))+
#     theme(legend.key.size=unit(.3,"cm"),
#           legend.box.background = element_rect(colour = "grey")) +
#     labs(x = 'Normalized Enrichment score')
#   
#   ggsave(paste0(path_temp, Tissue_t,"_GO_GENE_b27.pdf"), plot = g_temp, width = 12, height = 16, dpi = 300)
#   
# }

Tissues = t_dat$Tissue %>% unique()

b = t_dat %>% 
  filter(Tissue == "Colon_Transverse") %>%
  filter(!A == 0) %>%
  filter(!G == 0) %>%
  dplyr::mutate(coef_sign = ifelse(Coef > 0, 1, -1)) %>%
  dplyr::mutate(A_sign = A*coef_sign)

genelist = b$Coef

names(genelist) = b$gene

genelist = sort(genelist, decreasing = TRUE)

gs_A <-  gseGO(geneList=genelist, 
               ont ="BP", 
               keyType = "ENSEMBL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 1, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "fdr")

gs_res_A = data.frame(gs_A)

g = ggplot(gs_res_A)

g_temp = g +  geom_point(aes(x=NES,
                             y=reorder(Description,NES),
                             fill=p.adjust,
                             size=setSize),
                         shape=24,
                         stroke=.2)+
  theme_bw(base_size=16)+
  scale_fill_gradientn(colors=pal)+
  scale_y_discrete("")+
  scale_size(range=c(2,8))+
  theme(legend.key.size=unit(.3,"cm"),
        legend.box.background = element_rect(colour = "grey")) +
  labs(x = 'Normalized Enrichment score')
g_temp
  
gs_res = data.frame()

for (Tissue_t in Tissues) {
  print(Tissue_t)
  b = t_dat %>% 
    filter(Tissue == Tissue_t) %>%
    filter(!A == 0) %>%
    filter(!G == 0) %>%
    dplyr::mutate(coef_sign = ifelse(Coef > 0, 1, -1)) %>%
    dplyr::mutate(A_sign = A*coef_sign)
  
  genelist = b$A
  
  names(genelist) = b$gene
  
  genelist = sort(genelist, decreasing = TRUE)
  
  gs_A <-  gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 1, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "fdr")
  
  gs_res_A = data.frame(gs_A)
  
  genelist = b$G
  
  names(genelist) = b$gene
  
  genelist = sort(genelist, decreasing = TRUE)
  
  gs_G <-  gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 1, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "fdr")
  
  gs_res_G = data.frame(gs_G)
  
  gs_res_A$attribute = 'Age'
  
  gs_res_G$attribute = 'Genetics'
  
  gs_res_temp = rbind(gs_res_G, gs_res_A)
  
  gs_res_temp$Tissue_t = Tissue_t
  
  gs_res = rbind(gs_res, gs_res_temp)
}

write.csv2(gs_res, paste0(path_temp, 'gs_res_new.csv')) 

gs_res = vroom(paste0(path_temp, 'gs_res_new.csv'))

gs_res$pvalue

gs_res$pvalue = as.numeric(gsub(",", ".", gsub("\\.", "", gs_res$pvalue)))

gs_res$NES = as.numeric(gsub(",", ".", gsub("\\.", "", gs_res$NES)))

gs_res$p.adjust = as.numeric(gsub(",", ".", gsub("\\.", "", gs_res$p.adjust)))

gs_res_simp = gs_res %>% 
  dplyr::select(pvalue, attribute, 'Tissue_t', Description) %>%
  group_by(attribute, Tissue_t) %>%
  arrange(pvalue) %>%
  dplyr::mutate(exp_P=row_number()/n()) %>%
  ungroup() %>%
  dplyr::select(pvalue,exp_P,Description,attribute, Tissue_t) %>%
  dplyr::mutate(logP=-log(pvalue))

g = ggplot(gs_res_simp %>%
             filter(Tissue_t %in%
                      c("Whole_Blood", "Colon_Transverse",
                        "Liver", "Pancreas")) %>%
             mutate(Tissue = Tissue_t,
                    Tissue.type = Tissue %>% 
                      str_split("_") %>% 
                      sapply("[[", 1) %>% 
                      clean.subtypes,
                    Tissue.subtype = Tissue %>% 
                      str_split(., pattern = "_", n=2) %>% 
                      sapply(function(x){ifelse(length(x)>1, x[[2]], '')}) %>% 
                      clean.subtypes,
                    Tissue.neat = str_c(Tissue.type, " ", Tissue.subtype) %>% 
                      str_trim()
             ))

t_test = gs_res_simp %>% 
  dplyr::group_by(Tissue_t) %>%
  do(w=wilcox.test(logP~attribute,data=.,paired=FALSE)) %>%
  dplyr::summarize(Tissue_t, W=w$p.value, stat=w$statistic[1] )


g2 = g+geom_point(aes(x=-log(exp_P,10),y=logP,color=attribute),size=1)+
  theme_bw(base_size = 14)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=8))+
  geom_abline(aes(x = x), data = data.frame(x =c(0, 4)), slope = 1, intercept = 0) +
  labs(x = '-log(Expected)', y='-log(P values)', col = '') + 
  scale_color_brewer(palette="Set2") +
  facet_wrap(~Tissue.neat)

g2

ggsave(paste0(path_temp, "facet_pathway_qqplot_4.png"), plot = g2, width = 16, height = 12, dpi = 300, scale = 0.8)

Tissues_c7 = c('Esophagus_Mucosa', 'Colon_Transverse', 'Whole_Blood', 'Heart_Left_Ventricle',
                             'Esophagus_Gastroesophageal_Junction', 'Heart_Atrial_Appendage', 'Skin_Sun_Exposed_Lower_leg')
pal=pnw_palette("Bay",100)
gs_res %>% filter(Tissue_t == Tissue_t & attribute == "Age")%>% top_n(-20, wt = p.adjust)

for (Tissue_t in Tissues_c7) {
  

  g = ggplot(gs_res %>% 
               filter(Tissue_t == Tissue_t & attribute == "Age" & ONTOLOGY == 'BP') %>%
               top_n(-20, wt = p.adjust))

  g_temp = g +  geom_point(aes(x=NES,
                        y=reorder(Description,NES),
                        fill=p.adjust,
                        size=setSize),
                    shape=24,
                    stroke=.2)+
      theme_bw(base_size=16)+
      scale_fill_gradientn(colors=pal)+
      scale_y_discrete("")+
      scale_size(range=c(2,8))+
      theme(legend.key.size=unit(.3,"cm"),
            legend.box.background = element_rect(colour = "grey")) +
      labs(x = 'Normalized Enrichment score')
  g_temp
  
  ggsave(paste0(path_temp, Tissue_t,"_GO_AGE_b27.pdf"), plot = g_temp, width = 12, height = 16, dpi = 300)
}

a = gs_res_A %>% filter(qvalues < 0.05)
a
