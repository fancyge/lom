##############################################################################
### calculate icc based on shannon index
############################################################################## 
library(dplyr)
library(psych)

get_twin_icc <- function(data, tp, value_col) {
  
  df <- data %>%
    filter(Timepoint == tp) %>%
    mutate(
      pair_id = substr(Twin_ID, 1, 5),
      Sex_num = ifelse(Sex %in% c("M","Male","1"), 1, 0)
    )
  
  formula_str <- paste0(value_col, " ~ Sex_num")
  lm_model <- lm(as.formula(formula_str), data = df, na.action = na.exclude)
  df$adj_val <- residuals(lm_model)

  pair_counts <- df %>%
    count(pair_id) %>%
    filter(n == 2)
  
  df_wide <- df %>%
    semi_join(pair_counts, by = "pair_id") %>%
    group_by(pair_id) %>%
    arrange(Twin_ID) %>%
    summarise(
      val1 = first(adj_val),
      val2 = last(adj_val), 
      zygosity = first(zygosity),
      .groups = "drop"
    )
  
  mz_data <- df_wide %>% 
    filter(zygosity %in% c("MZ","MZ/MZ/MZ")) %>% 
    select(val1, val2)
  
  dz_data <- df_wide %>% 
    filter(zygosity %in% c("DZ","OSDZ","DZ/OSDZ/DZ")) %>% 
    select(val1, val2)
  
  icc_mz <- psych::ICC(mz_data)$results %>% filter(type == "ICC1")
  icc_dz <- psych::ICC(dz_data)$results %>% filter(type == "ICC1")
  
  res_table <- data.frame(
    Timepoint = tp,
    Zygosity  = c("MZ", "DZ"),
    N_pairs   = c(nrow(mz_data), nrow(dz_data)),
    ICC       = round(c(icc_mz$ICC, icc_dz$ICC), 3),
    CI_lower  = round(c(icc_mz$`lower bound`, icc_dz$`lower bound`), 3),
    CI_upper  = round(c(icc_mz$`upper bound`, icc_dz$`upper bound`), 3),
    p_value   = c(icc_mz$p, icc_dz$p)
  )
  
  return(res_table)
}

shannon <- read.csv("./inputs/shannon_map.csv")
icc_T1 <- get_twin_icc(shannon, "T1", "shannon_composition") #shannon_pwy_stratified_filtered; shannon_pwy_unstratified_noun
icc_T2 <- get_twin_icc(shannon, "T2", "shannon_composition")
icc_T3 <- get_twin_icc(shannon, "T3", "shannon_composition")

icc_all <- bind_rows(icc_T1, icc_T2, icc_T3)
print(icc_all)

