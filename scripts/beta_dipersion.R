####################################################################
#### test on beta dispersion across three timepoints
####################################################################
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr) 

## bacterial composition
abun_comp <- read.csv("./inputs/comp_metaphlan.csv")
abund <- abun_comp[, -c(1:4)]  
dist_bc_comp <- vegdist(abund, method = "bray") 
disp_bray <- betadisper(dist_bc_comp, meta$Timepoint)

disp_sp_tp <- data.frame(
   Timepoint = meta$Timepoint,
   Distance_to_Median = disp_bray $distances
 )

## species-stratified pwys
abun_func <- read.delim("./func_plsda_batch_corrected_t2b1_920_info.txt",header=T) 

func <- abun_func[,-c(1)] 
dist_mat <- dist(func, method = "euclidean") 
disp_func <- betadisper(dist_mat, meta$Timepoint) 

disp_fc_tp <- data.frame(
  Timepoint = meta$Timepoint,
  Distance_to_Median = disp_func $distances
)

## unifrac distance
unifrac_mat <- read.delim("./unifrac_merged_mpa3_profiles_920_log10.tsv") 
unifrac_dist <- as.dist(unifrac_mat) 
disp_unifrac <- betadisper(unifrac_dist, meta$Timepoint)

disp_unifrac_tp <- data.frame(
   Timepoint = meta$Timepoint,
   Distance_to_Median = disp_unifrac$distances
 )

disp_sp_tp$Metric <- "Bray-Curtis (Taxonomy)"
disp_fc_tp$Metric <- "Euclidean (Functions)"
disp_unifrac_tp$Metric <- "Weighted UniFrac"

combined_disp <- bind_rows(
    disp_sp_tp %>% select(Timepoint, Distance_to_Median, Metric),
    disp_fc_tp %>% select(Timepoint, Distance_to_Median, Metric),
    disp_unifrac_tp %>% select(Timepoint, Distance_to_Median, Metric) # Fixed the missing comma here
)

combined_disp$Metric <- factor(combined_disp$Metric, 
                               levels = c("Bray-Curtis (Taxonomy)", "Euclidean (Functions)","Weighted UniFrac"))

my_comparisons <- list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"))

ggplot(combined_disp, aes(x = Timepoint, y = Distance_to_Median, fill = Timepoint)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 0.8,width=0.6) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = c("T1" = "#424242", "T2" = "#E64626", "T3" = "#0148A4")) +
  theme_bw() +
  labs(
    y = "Distance to Group Median",
    x = "Timepoint"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", 
                     label = "p.signif",    
                     step.increase = 0.1,   
                     vjust = 0.8)                                
                     
                     
                     