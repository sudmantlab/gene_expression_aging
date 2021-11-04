library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(vroom)
library(dplyr)
library(tidyr)
library(stringr)


path = "Y:/analysis/Predixcan_eQTL/"

path_temp <- 'G:/Shared drives/sudmantlab/projects/agingGeneRegulation/figures/figure2x/fig2x_'

frame = vroom(paste0(path, "r2_b27.tsv"))

f = frame %>% select(Tissue) %>% unique()

f

tissues_c7 = c('Esophagus_Mucosa', 'Colon_Transverse', 'Whole_Blood', 'Heart_Left_Ventricle',
               'Esophagus_Gastroesophageal_Junction', 'Heart_Atrial_Appendage', 'Skin_Sun_Exposed_Lower_leg')

for (tissue in tissues_c7) {
  
}

df_genes = frame %>%
  filter(Tissue %in% tissues_c7) %>%
  rowwise() %>%
  mutate(gene = str_remove(gene, "\\..*")) %>%
  dplyr::select(Tissue, age, gene, R2) %>%
  unique() %>%
  spread(age, R2) %>%
  mutate(diff = young - old) %>%
  select(diff, gene, Tissue)

rownames(df_genes) = df_genes$gene

for(tissue in tissues_c7) {
  tissue
  gene_list_tis = df_genes %>% filter(Tissue == tissue)
  gene_list = gene_list_tis$diff
  names(gene_list) <- gene_list_tis$gene
  
  gene_list <- na.omit(gene_list) 
  
  gene_list = sort(gene_list, decreasing = TRUE)
  gs <-  gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  
  dotplot(gs, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle(tissue)
  
  ggsave(paste0(path_temp, tissue, '_GO.pdf'))
  
  d = data.frame(gs)
  
  
  
  write.csv2(d, paste0(path, tissue, '_GO.csv'))
  
  
}
library(PNWColors)

pal=pnw_palette("Bay",100)

for (tissue in tissues_c7) {
  d = read.csv2(paste0(path, tissue, '_GO.csv'))
  
  d_top = d %>% filter(p.adjust < 0.05 & enrichmentScore > 0.0) %>% top_n(10, -p.adjust)
  
  d_bot = d %>% filter(p.adjust < 0.05 & enrichmentScore < 0.0) %>% top_n(10, -p.adjust)
  
  g=ggplot(rbind(d_top, d_bot) %>% arrange(NES))
  
  g + geom_segment(aes(x=0,
                       xend=-enrichmentScore,
                       y=reorder(Description,enrichmentScore),
                       yend=reorder(Description,enrichmentScore)),
                   size=.3,
                   color="grey") +
    geom_point(aes(x=0,
                   y=reorder(Description,enrichmentScore)),
               size=.3,
               color="grey") +
    geom_point(aes(x=-enrichmentScore,
                   y=reorder(Description,enrichmentScore),
                   fill=pvalue,
                   size=setSize),
               shape=24,
               stroke=.2)+
    theme_bw(base_size=16)+
    scale_fill_gradientn(colors=pal)+
    scale_y_discrete("")+
    scale_size(range=c(2,8))+
    guides(size=FALSE)+
    theme(legend.position=c(0.82,0.85),
          legend.key.size=unit(.3,"cm"),
          legend.box.background = element_rect(colour = "grey")) +
    ggtitle(tissue) + 
    labs(x = 'Enrichment score')
  ggsave(paste0(path_temp, tissue, '_GO_nice.pdf'))
}

d = read.csv2(paste0(path, tissue, '_GO.csv'))
d$Tissue = tissue

for (tissue in tissues_c7) {
  temp = read.csv2(paste0(path, tissue, '_GO.csv'))
  temp$Tissue = tissue
  
  d = rbind(temp, d)
  
}

common = d %>% unique() %>%
  select(ID, Description, Tissue, NES) %>%
  spread(Tissue, NES) 

df = common[rowSums(is.na(common)) < 4, ]

ggplot(df) + geom_tile(y = Description, fill = NES)

na.omit(common)  
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
organism = org.Hs.eg.db


