get_pairwise_bray

##########################################################################################
####### get_pairwise_bray_sp
##########################################################################################

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr) 

data_abun <- read.csv("./inputs/abun_metaphlan_920.csv") # 734 species
rownames(data_abun) <- data_abun$Sample_ID 
abundance_matrix <- data_abun[, -c(1:7)] # Keep only the species columns

get_pairwise_bray <- function(tp, abun_mat, meta_df) {
  
  meta_tp <- meta_df %>% filter(Timepoint == tp)
  
  samples_to_keep <- intersect(rownames(abun_mat), meta_tp$Sample_ID)
  abun_tp <- abun_mat[samples_to_keep, ]
  
  bray_dist <- vegdist(abun_tp, method = "bray")
  
  dist_mat <- as.matrix(bray_dist)
  
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA 
  
  dist_df <- as.data.frame(as.table(dist_mat)) %>%
    rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
    filter(!is.na(Distance)) # Remove NAs (duplicates and self-comparisons)
  
  dist_df <- dist_df %>%
    left_join(meta_tp %>% select(Sample_ID, Family_ID, zygosity), 
              by = c("Sample1" = "Sample_ID")) %>%
    rename(Family1 = Family_ID, zygosity1 = zygosity)
  
  dist_df <- dist_df %>%
    left_join(meta_tp %>% select(Sample_ID, Family_ID, zygosity), 
              by = c("Sample2" = "Sample_ID")) %>%
    rename(Family2 = Family_ID, zygosity2 = zygosity)
  
  dist_df <- dist_df %>%
    mutate(
      Relationship = case_when(
        Family1 == Family2 & zygosity1 %in% c("MZ", "MZ/MZ/MZ") ~ "MZ Pair",
        Family1 == Family2 & zygosity1 %in% c("DZ", "OSDZ", "DZ/OSDZ/DZ") ~ "DZ Pair",
        Family1 != Family2 ~ "Unrelated Pair",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Relationship)) %>%
    mutate(Timepoint = tp)
  
  return(dist_df)
}

bray_T1 <- get_pairwise_bray("T1", abundance_matrix, meta)
bray_T2 <- get_pairwise_bray("T2", abundance_matrix, meta)
bray_T3 <- get_pairwise_bray("T3", abundance_matrix, meta)

bray_all <- bind_rows(bray_T1, bray_T2, bray_T3)

bray_all$Relationship <- factor(bray_all$Relationship, 
                                levels = c("MZ Pair", "DZ Pair", "Unrelated Pair"))
                                
                                
my_comparisons <- list(c("MZ Pair", "DZ Pair"), 
                       c("DZ Pair", "Unrelated Pair"), 
                       c("MZ Pair", "Unrelated Pair"))

summary_stats <- bray_all %>%
  group_by(Timepoint, Relationship) %>%
  summarise(
    Median_Distance = median(Distance),
    N_Pairs = n(),
    .groups = "drop"
  )
print(summary_stats)               


p <- ggplot(bray_all, aes(x = Relationship, y = Distance, fill = Relationship)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ Timepoint) +
  theme_bw() +
  scale_fill_manual(values = c("MZ Pair" = "#F8766D", 
                               "DZ Pair" = "#00BA38", 
                               "Unrelated Pair" = "#619CFF")) +
  labs(y = "Bray-Curtis Dissimilarity", x = "", 
       title = "Community Dissimilarity by Relationship") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     label = "p.signif") 

print(p)

                 