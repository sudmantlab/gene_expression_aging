# library(datasets)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
library(ggpubr)
library(tidyverse)
library(scales)
library(ggrepel)
library(ggsci)
library(dlookr)
# library(purrr)
library(vroom)
library(viridis)
# library(aplot)
library(cowplot)
# library(reshape2)
source("../Rscript/cleanSubType.R")

frame <- vroom('../analysis/JSD_All/All__JSD_result.csv.bz2')

path <- '../figures/figure_ryo/figures_adj_corr_no_sex/'
dir.create(path, recursive = T)
# frame

s = c('Cervix_Endocervix', 'Bladder', 'FallopianTube', 'Cervix_Ectocervix', 'Kidney_Cortex')

frame.clean <- frame %>%
  mutate(diff = old - young) %>%
  filter(!frame$Tissue %in% s) %>%
  mutate(
    sig = ifelse(`p-value` < 0.05,ifelse(diff > 0, "Significantly Increased", "Significantly Decreased"), "Not Significant"),
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
  "Significantly Increased" = "#D81B60", 
  "Significantly Decreased" = "#1E88E5", 
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
  # xlim(c(-1,1)) +
  # theme_void() + 
  theme_minimal() +
  theme(
    legend.position = "left",
    text = element_text(size=14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank()
  )

# g1  

g1.2 <- frame.clean %>% 
  ggplot(
    aes(
      x = diff, 
      y = reorder(Tissue.neat, diff), 
      col=sig
      )
         ) + 
  geom_point(size=5) + 
  scale_colour_manual(values = pal.sig, guide = guide_legend(title="")) + 
  labs(
    title = 'Variance difference among tissues',
    x = 'var(old) - var(young)'
  )  +
  theme_pubclean()+
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
   #scale_color_manual(values=c("#999999", "#56B4E9", "#E84826"))

g1.2

# # Using `aplot`
# g3 <- g2 %>% 
#   insert_left(g1)

# Using cowplot:
g1 <- plot_grid(g1.1+ylim2(g1.1), g1.2+ylim2(g1.1), align = "h", nrow = 1, ncol=2, axis="tb", rel_widths = c(1,2))
g1

ggsave(filename = paste0(path, "JSD_all.pdf"), plot = g3)