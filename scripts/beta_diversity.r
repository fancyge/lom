####################################################################################
##### beta diversity of timepoint groups. 
##### Plot with pcoa; statistics with PERMANOVA, anova, permutest
####################################################################################

library(vegan)       
library(ggplot2)    
library(dplyr)       
library(tibble)     
library(ape)
library(stringr)

############################## Bacterial composition ######################################## 
abun_comp <- read.csv("./inputs/comp_metaphlan.csv")

meta <- data.frame(Timepoint = abun_comp$Timepoint, Sex = abun_comp$Sex,Family_ID =str_extract(abun_comp$Sample_ID, "^T\\d+"))
abund <- abun_comp[, -c(1:4)]  

# calculate Bray-Curtis distance
dist_bc_comp <- vegdist(abund, method = "bray") 
pcoa_res_comp <- pcoa(dist_bc_comp)

### plot
pcoa_df_comp <- data.frame(
  Sample = 1:nrow(abun_comp),
  Timepoint = meta$Timepoint,
  Axis1 = pcoa_res_comp $vectors[,1],
  Axis2 = pcoa_res_comp $vectors[,2]
)

var_exp <- round(pcoa_res_comp$values$Relative_eig[1:2] * 100, 1) # Percent variance explained

ggplot(pcoa_df_comp, aes(x = Axis1, y = Axis2, color = Timepoint)) +
  geom_point(size = 0.5, alpha = 1) +
  stat_ellipse(aes(color = Timepoint), level = 0.95, show.legend = FALSE,linewidth = 0.2) +
  scale_colour_manual(values = c("T1" = "#424242", "T2" = "#E64626", "T3" = "#0148A4")) +
  labs(
    title = "Bacterial Composition", #PCoA (Bray-Curtis)
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Timepoints"
  ) +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    panel.background = element_blank(),
	panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = ""
  )
  
### stats

# PERMANOVA for significance, accounts for the repeated measures/families
permanova_comp <- adonis2(dist_bc_comp ~ Timepoint+Sex, data = meta, strata = meta$Family_ID, permutations = 999)

# multivariate dispersion 
disp_comp <- betadisper(dist_bc_comp, meta$Timepoint) 
anova(disp_comp)

# Permutation test for homogeneity of multivariate dispersions
perm_design <- how(nperm = 999, blocks = meta$Family_ID) 
permutest(disp_comp, permutations = perm_design) 

############################## Functional Diversity using filtered stratified pathways ################################# 
abun_func <- read.delim("./inputs/func_plsda.txt",header=T) # 2196, batch-corrected

func <- abun_func[,-c(1)] # data frame only, no mapping

# Euclidean distance for CLR data
dist_mat <- dist(func, method = "euclidean") 
pcoa_res_func <- pcoa(dist_mat)

### plot
plot_df <- data.frame(
  SampleID = rownames(func),
  Timepoint = meta$Timepoint,
  Axis1 = pcoa_res_func $vectors[,1],
  Axis2 = pcoa_res_func $vectors[,2]
)

ggplot(plot_df, aes(x = Axis1, y = Axis2, color = Timepoint)) +
  geom_point(size = 0.5, alpha = 1) +
  labs(
    title = "Functional Pathways", 
    x = paste0("PCoA1 (", round(pcoa_res_func $values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_res_func $values$Relative_eig[2] * 100, 1), "%)")
  ) +
  theme_minimal() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  )+
  theme(
    text = element_text(size = 10),
    legend.position = ""
  )+stat_ellipse(aes(color = Timepoint), level = 0.95, show.legend = FALSE,linewidth =0.2)+scale_colour_manual(values=c("T1"="#424242", "T2"="#E64626", "T3"="#0148A4"))+labs(color="Timepoints")

### stats 
permanova_func <- adonis2(dist_mat ~ Timepoint+Sex, data = meta, permutations = 999,strata = meta$Family_ID)

disp_func <- betadisper(dist_mat, meta$Timepoint) 
anova(disp_func)

perm_design <- how(nperm = 999, blocks = meta$Family_ID) 
permutest(disp_func, permutations = perm_design) 


## unifrac matrix
unifrac_mat <- read.delim("./inputs/unifrac_matrix.tsv") # distance data
unifrac_dist <- as.dist(unifrac_mat) 

permanova_unifrac <- adonis2(unifrac_dist ~ Timepoint+Sex, data = meta, permutations = 999,strata = meta$Family_ID)

disp_unifrac <- betadisper(unifrac_dist, meta$Timepoint) 
anova(disp_unifrac)

perm_design <- how(nperm = 999, blocks = meta$Family_ID) 
permutest(disp_unifrac, permutations = perm_design) 




