##################################################################################################
## Pairwise phylogenetic distances of S. mitis populations across developmental stages and relationship types
##################################################################################################
library(ape)
library(seqinr)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpubr)

tree <- read.tree("RAxML_bestTree.s__Streptococcus_mitis.StrainPhlAn3.tre")
dist_matrix <- cophenetic.phylo(tree)

dist_df <- as.data.frame(as.table(dist_matrix))

dist_df$Var1 <- as.character(dist_df$Var1)
dist_df$Var2 <- as.character(dist_df$Var2)

dist_clean <- dist_df %>%
  filter(Var1 < Var2) %>%
  rename(Distance = Freq)

dist_annotated <- dist_clean %>%
  left_join(meta %>% select(Sample_ID, Twin_ID, Family_ID, Timepoint), by = c("Var1" = "Sample_ID")) %>%
  rename(Twin_ID1 = Twin_ID, Family_ID1 = Family_ID, Time1 = Timepoint) %>%
  left_join(meta %>% select(Sample_ID, Twin_ID, Family_ID, Timepoint), by = c("Var2" = "Sample_ID")) %>%
  rename(Twin_ID2 = Twin_ID, Family_ID2 = Family_ID, Time2 = Timepoint) %>%
  mutate(
    Comparison_Group = case_when(
      Twin_ID1 != Twin_ID2 & Family_ID1 == Family_ID2 & Time1 == Time2 & Time1 == "T1" ~ "Twin (T1)",
      Twin_ID1 != Twin_ID2 & Family_ID1 == Family_ID2 & Time1 == Time2 & Time1 == "T2" ~ "Twin (T2)",
      Twin_ID1 != Twin_ID2 & Family_ID1 == Family_ID2 & Time1 == Time2 & Time1 == "T3" ~ "Twin (T3)",
      
      Twin_ID1 == Twin_ID2 & ((Time1 == "T1" & Time2 == "T2") | (Time1 == "T2" & Time2 == "T1")) ~ "Longitudinal (T1-T2)",
      Twin_ID1 == Twin_ID2 & ((Time1 == "T2" & Time2 == "T3") | (Time1 == "T3" & Time2 == "T2")) ~ "Longitudinal (T2-T3)",
      Twin_ID1 == Twin_ID2 & ((Time1 == "T1" & Time2 == "T3") | (Time1 == "T3" & Time2 == "T1")) ~ "Longitudinal (T1-T3)",
      
      TRUE ~ "Unrelated"
    )
  )

dist_annotated$Comparison_Group <- factor(
  dist_annotated$Comparison_Group, 
  levels = c("Unrelated", "Twin (T1)", "Twin (T2)", "Twin (T3)", 
             "Longitudinal (T1-T2)", "Longitudinal (T2-T3)", "Longitudinal (T1-T3)")
)

groups_to_test <- c("Twin (T1)", "Twin (T2)", "Twin (T3)", 
                    "Longitudinal (T1-T2)", "Longitudinal (T2-T3)", "Longitudinal (T1-T3)")

for (g in groups_to_test) {
  subset_data <- dist_annotated %>% filter(Comparison_Group %in% c(g, "Unrelated"))
  if(nrow(subset_data) > 0) {
    res <- wilcox.test(Distance ~ Comparison_Group, data = subset_data)
    cat(paste0("\n--- Wilcoxon Test: ", g, " vs Unrelated ---\n"))
    print(res)
  }
}

comparisons_list <- list(
  c("Unrelated", "Twin (T1)"),
  c("Unrelated", "Twin (T2)"),
  c("Unrelated", "Twin (T3)"),
  c("Unrelated", "Longitudinal (T1-T2)"),
  c("Unrelated", "Longitudinal (T2-T3)"),
  c("Unrelated", "Longitudinal (T1-T3)")
)

p <- ggplot(dist_annotated, aes(x = Comparison_Group, y = Distance, fill = Comparison_Group)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size = 0.5, width=0.6) +
  theme_minimal(base_size = 12) +
  labs(y = "Pairwise Phylogenetic Distance (Cophenetic)", x = "") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    legend.position = "none", 
    legend.title = element_blank(),
    
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text = element_text(color="black", size=9)
  ) +
  stat_compare_means(
    comparisons = comparisons_list, 
    method = "wilcox.test", 
    label = "p.signif" # Displays stars (*, **, ***) or raw p-values ("p.format")
  )





