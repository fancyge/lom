##############################################################################
## heatmap and violin plot of polymorphisc rates per species per tp 
##############################################################################

library(pheatmap)    
library(ComplexHeatmap)  
library(ggplot2)
library(grid)
library(gridExtra)
library(dendextend)
library(tibble)
library(tidyverse)
library(dplyr)
library(wesanderson)
library(viridis)

dp<-read.delim("~/Documents/PROJECTS/0_NIH_projects/1_within_species/PolymorphicRates/Rates_tables/922/all44.txt",sep="")

## heatmap using median
dp_sum <- dp %>% group_by(species, Time_point) %>%
  summarize(med_rate = median(polymorphic_rates_median_by_marker, na.rm = TRUE),
            .groups = "drop")

rate_mat <- dp_sum %>%
  pivot_wider(names_from = Time_point, values_from = med_rate) %>%
  column_to_rownames("species")

timepoint_distributions <- dp_sum %>%
  group_by(Time_point) %>%
  summarize(
    Q1 = quantile(med_rate, 0.25, na.rm = TRUE),
    Median = median(med_rate, na.rm = TRUE),
    Q3 = quantile(med_rate, 0.75, na.rm = TRUE),
    IQR = IQR(med_rate, na.rm = TRUE)
  )

row_dist <- dist(rate_mat, method = "euclidean")
hc       <- hclust(row_dist, method = "complete")
dend     <- as.dendrogram(hc)

dend2           <- dend            
dend2[[2]]      <- rev(dend2[[2]]) 

hc2 <- as.hclust(dend2)

heat_grob <- grid.grabExpr({
	pheatmap::pheatmap(
  	rate_mat,
 	color = viridis::viridis(100, option = "magma"),
	cluster_rows = hc2,
  	cluster_cols = FALSE,
  	fontsize_row  = 6,
  	fontsize_col  = 7,
    angle_col     = 0, 
    fontface = "italic"
)
})

## violin plot 
species_order2 <- hc2$labels[hc2$order]
dp2 <- dp %>% 
  mutate(
    species    = factor(species, levels = species_order2),
    Time_point = factor(Time_point, levels = c("T1","T2","T3"))
  )

p_violin <- ggplot(dp2, aes(Time_point, polymorphic_rates_median_by_marker, fill = Time_point)) +
  geom_violin(trim = FALSE, alpha = 0.7,size=0.2) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.size = 0.1) +
  scale_fill_manual(values = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")) +
  facet_wrap(~ species, ncol = 6, scales = "fixed") + 
  theme_bw(base_size = 6) +
  theme(
    axis.text.y    = element_text(size=6),
	axis.title.y     = element_blank(),
    axis.text.x      = element_text(size = 6),
    strip.text.x     = element_text(size = 7),
    panel.grid       = element_blank(),
    legend.position  = ""
  ) +
  labs(x = NULL, y = "Polymorphic rates")

violin_grob <- ggplotGrob(p_violin)




