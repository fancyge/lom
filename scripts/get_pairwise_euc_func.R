########################################################################
## get get_pairwise_euclidean distances from functional profiling 
## input is plsda correct func 
########################################################################

get_pairwise_euclidean <- function(tp, data_matrix, meta_df) {
  
  meta_tp <- meta_df %>% filter(Timepoint == tp)
  
  samples_to_keep <- intersect(rownames(data_matrix), meta_tp$Sample_ID)
  data_tp <- data_matrix[samples_to_keep, ]
  
  euc_dist <- dist(data_tp, method = "euclidean")
  
  dist_mat <- as.matrix(euc_dist)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA 
  
  dist_df <- as.data.frame(as.table(dist_mat)) %>%
    rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
    filter(!is.na(Distance))
  
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



