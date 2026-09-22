##############################################################################
#### run ace for individual features
############################################################################## 

library(dplyr)
library(tidyr)
library(stringr)
library(umx)

inv_norm_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

run_ace_feature <- function(dat, feature, tp, apply_int = TRUE) {
  message(sprintf("[%s] Running ACE for: %s", tp, feature))
  
  empty_res <- data.frame(
    Timepoint = tp,
    feature   = feature,
    A         = NA_real_,
    C         = NA_real_,
    E         = NA_real_,
    n_pairs   = 0L,
    p_value_A = NA_real_,
    p_value_C = NA_real_,
    p_value_E = NA_real_,
    stringsAsFactors = FALSE
  )
  
  df1 <- dat %>%
    filter(Timepoint == tp) %>%
    mutate(
      family_id = str_extract(Sample_ID, "^T\\d+"),
      twin_id   = str_extract(Sample_ID, "^T\\d+([A-Z])") %>% str_sub(-1),
      zygosity  = case_when(
        zygosity %in% c("MZ", "MZ/MZ/MZ") ~ "MZ",
        zygosity %in% c("DZ", "OSDZ", "DZ/OSDZ/DZ") ~ "DZ",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(twin_id %in% c("A", "B"), !is.na(zygosity)) %>%
    distinct(family_id, twin_id, .keep_all = TRUE)
  
  raw_vals <- df1[[feature]]
  
  if (apply_int) {
    df1$value <- as.numeric(scale(inv_norm_transform(raw_vals)))
  } else {
    df1$value <- as.numeric(scale(raw_vals))
  }
  
  df_wide <- df1 %>%
    pivot_wider(
      id_cols     = c(family_id, zygosity),
      names_from  = twin_id,
      values_from = c(value, Sex)
    ) %>%
    filter(!is.na(value_A), !is.na(value_B))
  
  n_pairs <- nrow(df_wide)
  empty_res$n_pairs <- n_pairs
  
  if (n_pairs < 10) return(empty_res)
  
  df_wide <- df_wide %>%
    rename(
      alpha1 = value_A,
      alpha2 = value_B,
      sex1   = Sex_A,
      sex2   = Sex_B
    ) %>%
    mutate(
      sex1 = ifelse(sex1 %in% c("M", "Male", "1"), 1, 0),
      sex2 = ifelse(sex2 %in% c("M", "Male", "1"), 1, 0)
    )
  
  m1 <- tryCatch(
    umxACE(
      selDVs  = "alpha",
      selCovs = "sex",
      mzData  = df_wide %>% filter(zygosity == "MZ"),
      dzData  = df_wide %>% filter(zygosity == "DZ"),
      sep     = "",
      autoRun = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(m1)) return(empty_res)
  
  su <- tryCatch(umxSummary(m1, std = TRUE), error = function(e) NULL)
  if (is.null(su)) return(empty_res)
  
  pct_A <- round(su[, "a1"]^2 * 100, 1)
  pct_C <- round(su[, "c1"]^2 * 100, 1)
  pct_E <- round(su[, "e1"]^2 * 100, 1)
  
  get_nested_p <- function(sub_model) {
    if (is.null(sub_model)) return(NA_real_)
    tryCatch({
      comp <- umxCompare(m1, sub_model)
      comp$p[2]
    }, error = function(e) NA_real_)
  }
  
  ceModel <- tryCatch(umxModify(m1, update = "a_r1c1", free = FALSE), error = function(e) NULL)
  aeModel <- tryCatch(umxModify(m1, update = "c_r1c1", free = FALSE), error = function(e) NULL)
  eModel  <- tryCatch(umxModify(m1, update = c("a_r1c1", "c_r1c1"), free = FALSE), error = function(e) NULL)
  
  data.frame(
    Timepoint = tp,
    feature   = feature,
    A         = pct_A,
    C         = pct_C,
    E         = pct_E,
    n_pairs   = n_pairs,
    p_value_A = get_nested_p(ceModel),
    p_value_C = get_nested_p(aeModel),
    p_value_E = get_nested_p(eModel),
    stringsAsFactors = FALSE
  )
}