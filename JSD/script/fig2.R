library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(scales)
library(ggrepel)
library(ggsci)
library(purrr)
library(vroom)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(broom)
library(plyr)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(see)
library(gghalves)
library(aplot)
library(stringr)
library(latex2exp)
library(DEGreport)

setwd("Z:/sudmant_lab/gene_expression_variance/Rscript")
source('../Rscript/cleanSubType.R')
source('../Rscript/split_violin.R')


### Load dataframes and set up ###

p_cutoff <- 0.00001

remove <- c("Brain_Amygdala", "Brain_Substantianigra", "Brain_Spinalcordcervicalc_1", "Kidney_Cortex",
            "Bladder", "Cervix_Endocervix", "FallopianTube", "Cervix_Ectocervix"
)

one = all_dist_frame %>% filter(value == 0) %>% select(Tissue, age) %>% unique %>% filter(Tissue %in% gt50)

frame <- vroom('Y:/analysis/JSD_All/All_JSD_result.csv') %>% 
  filter(!Tissue %in% remove)

all_dist_frame <- vroom('../analysis/JSD_All/All__dist.csv.bz2')%>% 
  filter(!Tissue %in% remove)

dist_frame <- vroom('Y:/analysis/JSD_binned/All_dist.csv')%>% 
  filter(!Tissue %in% remove)

path_temp <- 'G:/Shared drives/sudmantlab/projects/agingGeneRegulation/figures/figure2/fig2_'
#path_temp <- '../figures/figures_adj_corr/JSD/'
dir.create(path_temp, recursive = T)

path_fig <- '..figures/figure2/fig2_'

path <- path_temp
dir.create(path, recursive = T)

tissue_names <-  vroom('../analysis/misc/selected_tissue_abbreviations2.tsv')

tissue_names

frame$Tissue %>% unique %>% sort

frame$Tissue_ab <- mapvalues(frame$Tissue, tissue_names$SMTSD, tissue_names$tissue)
dist_frame$Tissue_ab <- mapvalues(dist_frame$Tissue, tissue_names$SMTSD, tissue_names$tissue)
all_dist_frame$Tissue_ab <- mapvalues(all_dist_frame$Tissue, tissue_names$SMTSD, tissue_names$tissue)

stats_frame <- dist_frame %>% 
  group_by(Tissue) %>%
  do(tidy(lm(value ~ age, data = .)))


stats_frame.clean <- stats_frame %>%
  filter(term == 'age') %>%
  filter(!Tissue %in% remove) %>%
  mutate(
    sig = ifelse(`p.value` < p_cutoff,
                 ifelse(estimate > 0, 
                        "Sig. Increased", 
                        "Sig. Decreased"),
                 "Not Significant"),
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

frame.clean <- frame %>%
  mutate(diff = old - young) %>%
  mutate(
    sig = ifelse(`p-value` < p_cutoff, ifelse(diff > 0, "Sig. Increased", "Sig. Decreased"), "Not Significant"),
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



### Plot Correlation ###

joint.frame <- left_join(stats_frame.clean %>% dplyr::rename(sig.estimate=sig),
                         frame.clean %>% dplyr::select(Tissue, diff, sig.diff = sig), 
                         by="Tissue")

#write.csv(joint.frame %>% select(Tissue, estimate, diff) %>% arrange(-estimate), paste0(path, 'B_summary.csv'))



gt50 = c('WholeBlood','Stomach', 'Colon_Sigmoid', 'Esophagus_GastroesophagealJunction', 'Colon_Transverse',
            'Artery_Aorta',
            'Heart_AtrialAppendage',
            'Breast_MammaryTissue',
            'Prostate',
            'Heart_LeftVentricle',
            'Cells_Culturedfibroblasts',
            'Esophagus_Muscularis',
            'Adipose_VisceralOmentum',
            'Esophagus_Mucosa',
            'Lung',
            'Skin_NotSunExposedSuprapubic',
            'Nerve_Tibial',
            'Thyroid',
            'Adipose_Subcutaneous',
            'Artery_Tibial',
            'Testis',
            'Skin_SunExposedLowerleg',
            'Muscle_Skeletal',
            'Pancreas',
            'Liver',
            'Artery_Coronary',
            'AdrenalGland'
         )

gt50

joint.frame$Tissue_ab <- mapvalues(joint.frame$Tissue, tissue_names$SMTSD, tissue_names$tissue)

joint.frame = joint.frame %>% filter(Tissue %in% gt50)

save(joint.frame, file="joint_frame.Rda")



write.csv2(joint.frame, "../analysis/JSD_All/JSD_joint_result.csv")

R <- cor.test(joint.frame$estimate, joint.frame$diff, use='all.obs', method="spearman")

R$estimate

R

R <- signif(R$estimate, digits=4)

R

tissues <- joint.frame %>% filter(Tissue %in% gt50)%>% arrange(-estimate) %>% pull(Tissue.neat) %>% unique
# pal.tissue <- pal_d3(palette = "category20")(length(tissues)) %>% 
#   set_names(., tissues)

saveRDS(pal.tissue, file = '../analysis/pal_tissue.Rdata')

a <- readRDS('../analysis/pal_tissue.Rdata')

a

pal.tissue = a

g5.main <- joint.frame %>% filter(Tissue %in% gt50) %>%
  ggplot(
    aes(
      x=estimate,
      y=diff,
      label=Tissue_ab
    )
  ) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +  
  geom_text_repel(force = 1, nudge_y = 0.001, size=5) + 
  geom_point(aes(fill=Tissue.neat), size=5, pch=21, colour='black') +
  labs(
    x="JSD slope",
    y=TeX('$\\Delta{JSD}$ (old - young )')
  )+
  theme_bw(base_size=16)+
  theme(
    legend.position = "None",
    legend.box = "vertical",
    text = element_text(size=14),
    axis.title = element_text(size=16),
    axis.text = element_text(size=16),
    plot.margin = unit(c(0, 1, 0, 0), "cm")
  ) + 
  scale_fill_manual(values = pal.tissue
                    ) + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  #geom_errorbar(aes(xmin=estimate-std.error, xmax=estimate+std.error)) +
  annotate(geom='text', x=-0.4e-3, y=0.015, label=paste("R = ", R), size=5 )

# g5.density <- ggMarginal(g5.main, type="density")
# g5.densigram <- ggMarginal(g5.main, type="densigram", yparams = list(binwidth = 0.003), xparams = list(binwidth = 0.8e-4))
# 
# g5.densigram.gridlines <- (g5.main +
#   theme(
#     panel.grid.major = element_line(linetype = "dotted", color="grey")
#     # axis.line.y = element_line(linetype = "dotted", color="lightgrey")
#   )) %>% ggMarginal(type="densigram", yparams = list(binwidth = 0.003), xparams = list(binwidth = 0.8e-4))
# g5 <- g5.densigram.gridlines

g5.main

g5 = g5.main

ggsave(filename = paste0(path, "tissues_scatterPlot_slope-v-diff.pdf"), plot=g5, width = 14, height = 12, scale=0.5)

joint.frame %>%
  mutate(up = ifelse(estimate > 0 & diff > 0, 'up', 'down')) %>%
  filter(up == 'up') %>%
 summarize(es_median = median(estimate), di_median = median(diff))

joint.frame %>%
    mutate(up = ifelse(estimate > 0 & diff > 0, 'up', 'down')) %>%
    filter(up == 'down') %>%
    summarize(es_median = median(estimate), di_median = median(diff))
  
### Plot Raindrop ###

all_dist_frame.clean <- all_dist_frame %>% 
  # filter(Tissue %in% select_tissues) %>% 
  mutate(
    age = factor(age %>% str_to_title(), levels = c("young", "old") %>% str_to_title()),
    Tissue.type = Tissue %>% 
      str_split("_") %>% 
      sapply("[[", 1) %>% 
      clean.subtypes,
    Tissue.subtype = Tissue %>% 
      str_split(., pattern = "_", n=2) %>% 
      sapply(function(x){ifelse(length(x)>1, x[[2]], '')}) %>% 
      clean.subtypes%>% 
      str_wrap(width="17"),
    Tissue.neat = str_c(Tissue.type, "\n", Tissue.subtype) %>% 
      str_trim()
  )

all_dist_frame.clean.summary <- all_dist_frame.clean %>% 
  dplyr::group_by(age, Tissue.neat, Tissue.type) %>% 
  dplyr::summarize(
    value.mean = value %>% mean(na.rm=T),
    value.median = value %>% median(na.rm=T),
    value.sd = value %>% sd(na.rm=T),
    value.n = n(),
    value.ci = qt(0.95,df=value.n-1)*value.sd/sqrt(value.n),
    value.ci_low = value.mean - value.ci,
    value.ci_high = value.mean + value.ci
  ) %>% 
  ungroup()

all_dist_frame.ageChange <- all_dist_frame.clean.summary %>% 
  dplyr::select(Tissue.neat, age, value.mean) %>% 
  pivot_wider(id_cols=Tissue.neat, names_from=age, values_from=value.mean) %>% 
  mutate(dvar = Old-Young) %>% 
  arrange(-dvar)

levels.var.tissue <- all_dist_frame.ageChange %>% 
  dplyr::pull(Tissue.neat)

all_dist_frame.clean$Tissue.neat %>% factor(., levels=levels.var.tissue)
all_dist_frame.clean.summary$Tissue.neat %>% factor(., levels=levels.var.tissue)

# Select 8 tissues of interest
## Criteria: abs(dvar) > 0.01
select_tissues <- c(
  # Var Increases with Age
 # 'Brain\nCortex',                               # Young 0.544 --> 0.574 Old #1 
 # 'Brain\nFrontal Cortex\n(Brodmann A. 9)',                     # Young 0.567 --> 0.653 Old #2
  'WholeBlood', 
  'Heart_LeftVentricle',
  'Stomach',
  'Colon_Transverse',
  'Prostate'
  #'Adipose_VisceralOmentum'
  )

select_dist_frame.clean <- all_dist_frame.clean %>% 
  filter(
    Tissue %in% select_tissues)
    # Tissue.neat %in% pal.tissue.2
  # ) %>% 
  # mutate(
  #   Tissue.neat = droplevels(Tissue.neat)
  # ) %>% mutate(up_down = ifelse(''))

select_dist_frame.clean.summary <- all_dist_frame.clean.summary %>% 
  filter(
    Tissue.neat %in% (select_dist_frame.clean$Tissue.neat %>% unique)
  )
  # ) %>% 
  # mutate(
  #   Tissue.neat = droplevels(Tissue.neat)
  # )

pal.tissue

# select_dist_frame.clean %>% select(Tissue.neat) %>% unique()
# 
# pal.tissue["Skin\nLower Leg (Sun\nExposed)"] = "#DBDB8DFF"

g6 <- 
  select_dist_frame.clean %>% filter(!value == 0.0) %>%
  ggplot(
    aes(
      fill = Tissue.neat, 
      y = value,
      x = age,
    )
  ) + 
  geom_half_violin(
    aes(fill=Tissue.neat)
  ) +
  geom_point(
    aes(x=factor(age) %>% as.numeric() %>% add(0.20)),
    position = position_jitter(width = .15),
    size = .25,
    alpha=0.25
  )+
  geom_pointrange(
    data = select_dist_frame.clean.summary,
    mapping = aes(
      x=age,
      ymin=value.ci_low,
      ymax=value.ci_high,
      y=value.mean
    ),
    color="red",
    size = 0.25
  )+
  geom_ribbon(
    data = select_dist_frame.clean.summary,
    mapping = aes(
      x = as.numeric(age),
      ymin=value.ci_low,
      ymax=value.ci_high,
      y=value.mean
    ),
    color="red",
    size = 0.5,
    alpha=0.25,
    inherit.aes = F
  )+
  facet_wrap(
    ~Tissue.neat, 
    nrow=1,
    drop = T
  ) + 
  scale_fill_manual(
    values = pal.tissue
  )+
  labs(
    x="Age Class",
    y="JSD",
    fill = "Tissue\nType"
  ) +
  theme_bw(base_size=24)+
  theme(
    legend.position = "None",
    legend.box = "vertical",
    text = element_text(size=24),
    axis.title.y = element_text(angle=0),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 24),
    axis.title = element_text(size=24),
    axis.text = element_text(size=20),
    strip.background = element_rect(fill="white", color="white")
  )

g6

ggsave(paste(path, "dist_tissues_facet_8.png", sep = ""), plot=g6, width = 16, height = 9)

# g + 
#   geom_violin(aes(fill = Tissue_ab),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
#   geom_point(aes(x = as.numeric(1)-.15, y = value, colour = Tissue_ab),position = position_jitter(width = .05), size = 0.1, shape = 1,alpha=0.8)+
#   geom_boxplot(aes(x = 1, y = value, fill = Tissue_ab),outlier.shape = NA, alpha = .5, width = .1, colour = "black",size=.1) +
#   facet_wrap(Tissue_ab~age)


var.test(values ~ groups, data, 
         alternative = "two.sided")

### Variation difference ### 
     
all_dist_frame.ftest <- all_dist_frame %>% 
  mutate(
    age = factor(age %>% str_to_title(), levels = c("young", "old") %>% str_to_title()),
    Tissue.type = Tissue %>% 
      str_split("_") %>% 
      sapply("[[", 1) %>% 
      clean.subtypes,
    Tissue.subtype = Tissue %>% 
      str_split(., pattern = "_", n=2) %>% 
      sapply(function(x){ifelse(length(x)>1, x[[2]], '')}) %>% 
      clean.subtypes%>% 
      str_wrap(width="17"),
    Tissue.neat = str_c(Tissue.type, "\n", Tissue.subtype) %>% 
      str_trim()
  ) %>% 
  select(-`Unnamed: 0`, -Tissue, -metric, -Tissue_ab) %>%
  nest(data = c(-Tissue.neat, -Tissue.type, -Tissue.subtype)) %>% 
  mutate(
    test = map(data, ~ fligner.test(.$value, .$age)), # S3 list-col
    tidied = map(test, tidy, conf.int=T)
  ) %>% 
  unnest(tidied)


g8.1 <-
  all_dist_frame.ftest %>% 
  filter(!Tissue.type %in% setdiff(all_dist_frame.ftest$Tissue.type, pal.tissue.2 %>% names)) %>% 
  ggplot(
    aes(
      x = 0, 
      y = reorder(Tissue.neat, statistic), 
      # col=p.value<p_cutoff,
      fill=Tissue.type,
      label=Tissue.type,
      width = 0.25
    )
  ) +
  scale_fill_manual(
    values = pal.tissue.2, 
    guide=guide_legend(
      title = "Tissue", 
      title.theme = element_text(face = "bold"),
      title.position = "top", 
      vjust=0.5,
      ncol=1,
      override.aes = list(shape=22, color="white")
    )
  )+
  geom_tile(color="white") + 
  geom_text() +
  theme_minimal() +
  theme(
    # legend.position = "left",
    legend.position = "none",
    text = element_text(size=14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank()
  )

g8.1

g8.2 <- all_dist_frame.ftest %>% 
  filter(!Tissue.type %in% setdiff(all_dist_frame.ftest$Tissue.type, pal.tissue.2 %>% names)) %>% 
  ggplot(
    aes(
      x = statistic, 
      y = reorder(Tissue.neat, statistic), 
      col=p.value<=p_cutoff
    )
  ) + 
  geom_point(size=5) + 
  geom_vline(xintercept = -log10(p_cutoff), lty="dashed", size = 1.5)+
  scale_colour_startrek()+ 
  labs(
    #title = 'Significance of variation differences between young and old tissues between individuals',
    x = 'Fligner Killeen Test'
  )  +
  theme_bw()+
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    legend.box = "vertical",
    text = element_text(size=14),
    plot.title = element_text(size = 18),
    axis.title = element_text(size=20),
    axis.text.x = element_text(size=16),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

g8 <- plot_grid(g8.1+ylim2(g8.2), g8.2+ylim2(g8.2), align = "h", nrow = 1, ncol=2, axis="tb", rel_widths = c(1,2))
g8


ggsave(filename = paste0(path, "tissues_var.pdf"), plot=g8, width = 8, height = 9)


############################################################################################################

############### Possible supplement figure code below ######################################################

############################################################################################################
### Plot Rank Tissues by variance difference ###

frame.clean <- frame %>%
  mutate(diff = old - young) %>%
  mutate(
    sig = ifelse(`p-value` < p_cutoff, ifelse(diff > 0, "Sig. Increased", "Sig. Decreased"), "Not Significant"),
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

pal.sig <- c(
  "Sig. Increased" = "#D81B60", 
  "Sig. Decreased" = "#1E88E5", 
  "Not Significant" = "#000000"
)

tissues <- frame.clean %>% arrange(-diff) %>% pull(Tissue.type) %>% unique
# pal.tissue.1 <- c(
#   "#c3a5b4", "#5d4c86", "#fcff5d", "#7dfc00", "#0ec434", "#228c68", 
#   "#8ad8e8", "#235b54", "#29bdab", "#3998f5", "#37294f", "#277da7", 
#   "#3750db", "#f22020", "#991919", "#ffcba5", "#e68f66", "#c56133", 
#   "#96341c", "#632819", "#ffc413", "#f47a22", "#2f2aa0", "#b732cc", 
#   "#772b9d", "#f07cab", "#d30b94"#, "#201923", "#edeff3", "#946aa2", "#ffffff"
# ) %>% 
#   set_names(., tissues)

pal.tissue.2 <- colorRampPalette(pal_d3(palette = "category20")(20))(length(tissues)) %>% 
  set_names(., tissues)

# # For future use of pal.tissue.2:
# sapply(names(pal.tissue.2), function(x, p=pal.tissue.2){sprintf('"%s" = "%s"', x, p[[x]])}) %>% paste0(collapse=",\n") %>% sprintf("c(%s)", .) %>% cat

g1.1 <-
  frame.clean %>% 
  ggplot(
    aes(
      x = 0, 
      y = reorder(Tissue.neat, diff), 
      # col=Tissue.type,
      fill=Tissue.type,
      label=Tissue.subtype %>% str_remove("\\(.*\\)") %>% str_trunc(side="center", width=15),
      width = 0.25
    )
  ) +
  scale_fill_manual(
    values = pal.tissue.2, 
    guide=guide_legend(
      title = "Tissue", 
      title.theme = element_text(face = "bold"),
      title.position = "top", 
      vjust=0.5,
      ncol=1,
      override.aes = list(shape=22, color="white")
    )
  )+
  geom_tile(color="white") + 
  geom_text() + 
  theme_minimal() +
  theme(
    # legend.position = "left",
    legend.position = "none",
    text = element_text(size=14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank()
  )

g1.2 <- frame.clean %>% 
  ggplot(
    aes(
      x = diff, 
      y = reorder(Tissue.neat, diff), 
      col=sig
    )
  ) + 
  geom_point(size=5) + 
  geom_vline(xintercept = 0, lty="dashed", size = 1.5)+
  scale_colour_manual(values = pal.sig, guide = guide_legend(title="")) + 
  labs(
    title = 'Variance difference among tissues',
    x = 'var(old) - var(young)'
  )  +
  theme_bw()+
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    legend.box = "vertical",
    text = element_text(size=14),
    plot.title = element_text(size = 18),
    axis.title = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# g1.2

g1 <- plot_grid(g1.1+ylim2(g1.1), g1.2+ylim2(g1.1), align = "h", nrow = 1, ncol=2, axis="tb", rel_widths = c(1,2))
g1

ggsave(plot= g1, filename = paste0(path, "tissues_rank_var_diff.new.pdf"), width = 6, height = 8)

### Plot tissue rank by slope ###

linear_dist <- dist_frame %>% 
  group_by(Tissue) %>%
  do(model = lm(value ~ age, data = .))

stats_frame <- tidy(linear_dist, model,
                    conf.int = T,
                    conf.level = 0.95)

stats_frame.clean <- stats_frame %>%
  filter(term == 'age') %>%
  filter(!Tissue %in% remove) %>%
  mutate(
    sig = ifelse(`p.value` < p_cutoff,
                 ifelse(estimate > 0, 
                        "Sig. Increased", 
                        "Sig. Decreased"),
                 "Not Significant"),
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

g2.1 <-
  stats_frame.clean %>% 
  ggplot(
    aes(
      x = 0, 
      y = reorder(Tissue.neat, estimate), 
      # col=Tissue.type,
      fill=Tissue.type,
      label=Tissue.subtype %>% str_remove("\\(.*\\)") %>% str_trunc(side="center", width=15),
      width = 0.25
    )
  ) +
  scale_fill_manual(
    values = pal.tissue.2, 
    guide=guide_legend(
      title = "Tissue", 
      title.theme = element_text(face = "bold"),
      title.position = "top", 
      vjust=0.5,
      ncol=1,
      override.aes = list(shape=22, color="white")
    )
  )+
  geom_tile(color="white") + 
  geom_text() + 
  theme_minimal() +
  theme(
    legend.position = "left",
    text = element_text(size=14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank()
  )

g2.2 <- stats_frame.clean %>% 
  ggplot(
    aes(
      x = estimate, 
      xmin = conf.low,
      xmax = conf.high,
      y = reorder(Tissue.neat, estimate), 
      col=sig
    )
  ) + 
  geom_linerange(alpha=0.5, color="black") +
  geom_point(size=5) + 
  geom_vline(xintercept = 0, lty="dashed", size = 1.5)+
  scale_colour_manual(values = pal.sig, guide = guide_legend(title="")) + 
  labs(
    title = 'Variance difference among tissues',
    x = 'Slope'
  )  +
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    text = element_text(size=14),
    plot.title = element_text(size = 18),
    axis.title = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

g2 <- plot_grid(g2.1+ylim2(g2.1), g2.2+ylim2(g2.1), align = "h", nrow = 1, ncol=2, axis="tb", rel_widths = c(1,2))
g2

ggsave(filename = paste0(path, "tissues_rank_var_slope.pdf"), plot=g2, width = 10, height = 3.5)


g8v2.1 <-
  all_dist_frame.ftest %>% 
  filter(!Tissue.type %in% setdiff(all_dist_frame.ftest$Tissue.type, pal.tissue.2 %>% names)) %>% 
  ggplot(
    aes(
      x = 0, 
      y = reorder(frame.clean$Tissue.neat, frame.clean$diff), #reorder(Tissue.neat, statistic), 
      # col=p.value<p_cutoff,
      fill=Tissue.type,
      label=Tissue.subtype %>% str_remove("\\(.*\\)") %>% str_trim() %>% str_trunc(side="center", width=15),
      width = 0.25
    )
  ) +
  scale_fill_manual(
    values = pal.tissue.2, 
    guide=guide_legend(
      title = "Tissue", 
      title.theme = element_text(face = "bold"),
      title.position = "top", 
      vjust=0.5,
      ncol=1,
      override.aes = list(shape=22, color="white")
    )
  )+
  geom_tile(color="white") + 
  geom_text() +
  theme_minimal() +
  theme(
    # legend.position = "left",
    legend.position = "none",
    text = element_text(size=14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank()
  )

g8v2.2 <- all_dist_frame.ftest %>% 
  filter(!Tissue.type %in% setdiff(all_dist_frame.ftest$Tissue.type, pal.tissue.2 %>% names)) %>% 
  ggplot(
    aes(
      x = statistic, 
      y = reorder(frame.clean$Tissue.neat, frame.clean$diff), #reorder(Tissue.neat, statistic), 
      col=p.value<=p_cutoff
    )
  ) + 
  geom_point(size=5) + 
  geom_vline(xintercept = -log10(p_cutoff), lty="dashed", size = 1.5)+
  scale_colour_startrek()+ 
  labs(
    title = 'Significance of variation differences between young and old tissues between individuals',
    x = 'Fligner Killeen Test'
  )  +
  theme_bw()+
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    legend.box = "vertical",
    text = element_text(size=14),
    plot.title = element_text(size = 18),
    axis.title = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

g9v2 <- plot_grid(g1.2+theme(legend.position = "none") + ylim2(g1.2), 
                  g1.1+ylim2(g1.2), 
                  g8v2.1+ylim2(g8v2.2), 
                  g8v2.2+theme(legend.position = "none")+ylim2(g8v2.2), 
                  align="h", nrow=1, ncol=4, axis="tb", rel_widths = c(2,1,1,2)) + theme(legend.position = "bottom") 


g7.1 <-
  all_dist_frame.ftest %>% 
  filter(!Tissue.type %in% setdiff(all_dist_frame.ftest$Tissue.type, pal.tissue.2 %>% names)) %>% 
  ggplot(
    aes(
      x = 0, 
      y = reorder(Tissue.neat, -log10(p.value)), 
      # col=Tissue.type,
      fill=Tissue.type,
      # label=Tissue.subtype %>% str_remove("\\(.*\\)") %>% str_trunc(side="center", width=15),
      width = 0.25
    )
  ) +
  scale_fill_manual(
    values = pal.tissue.2, 
    guide=guide_legend(
      title = "Tissue", 
      title.theme = element_text(face = "bold"),
      title.position = "top", 
      vjust=0.5,
      ncol=1,
      override.aes = list(shape=22, color="white")
    )
  )+
  geom_tile(color="white") + 
  # geom_text() +
  theme_minimal() +
  theme(
    legend.position = "left",
    text = element_text(size=14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank()
  )

g7.1

g7.2 <- all_dist_frame.ftest %>% 
  filter(!Tissue.type %in% setdiff(all_dist_frame.ftest$Tissue.type, pal.tissue.2 %>% names)) %>% 
  ggplot(
    aes(
      x = -log10(p.value), 
      y = reorder(Tissue.neat, -log10(p.value)), 
      col=p.value<=p_cutoff
    )
  ) + 
  geom_point(size=5) + 
  geom_vline(xintercept = -log10(p_cutoff), lty="dashed", size = 1.5)+
  scale_colour_startrek()+ 
  labs(
    #title = 'Significance of variation differences between young and old tissues between individuals',
    x = '-log(P-value)'
  )  +
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    text = element_text(size=14),
    plot.title = element_text(size = 18),
    axis.title = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

g7 <- plot_grid(g7.1+ylim2(g7.2), g7.2+ylim2(g7.2), align = "h", nrow = 1, ncol=2, axis="tb", rel_widths = c(1,2))
g7

