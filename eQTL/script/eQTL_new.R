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

frame = vroom("Y:/analysis/JSD_eQTL/Esophagus_Gastroesophageal_Junction_lm_results_new.tsv")

names(frame)

t_sum = frame %>% 
  mutate(P=`p_val_genotype[T.0|1]`) %>%
  group_by(Tissue,age) %>%
  arrange(P) %>%
  dplyr::mutate(exp_P=row_number()/n(),
                mu_P = mean(P)) %>%
  ungroup() %>%
  dplyr::select(P,exp_P,gene,Tissue,age) %>%
  mutate(logP=-log(P)) %>%
  mutate(logexpP = -log(exp_P)) %>%
  arrange(age)

t_label = t_sum %>% 
  do(w=wilcox.test(logP~age,data=.,paired=FALSE)) %>%
  dplyr::summarize(W=w$p.value, stat=w$statistic[1] )

a = wilcox.test(logP~age,data=t_sum,paired=FALSE)

a$p.value

g = ggplot(t_sum)

g + geom_point(aes(x = logexpP, y = logP, color = age))+ 
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  geom_text(aes(x = 2.5, y = 90, label = "P = 7.0e-35"), size = 8) +
  ggtitle("Esophagus_Gastroesophageal_Junction")

