library(umx)
library(dplyr)

#########################################################################
# Function 1: Run ACE model for a specific timepoint (Covariate: Sex)
#########################################################################
run_ace_timepoint <- function(data, tp, value_col){

  df <- data %>%
    filter(Timepoint == tp) %>%
    select(Sample_ID, Twin_ID, all_of(value_col), zygosity, Sex, Family_ID) %>%
    mutate(pair_id = substr(Twin_ID, 1, 5))
  
  df[[value_col]] <- as.numeric(scale(df[[value_col]]))
  
  # 3. Keep only complete pairs
  pair_counts <- df %>%
    count(pair_id) %>%
    filter(n == 2)
  
  df_wide <- df %>%
    semi_join(pair_counts, by = "pair_id") %>%
    group_by(pair_id) %>%
    arrange(Twin_ID) %>%
    summarise(
      alpha1   = first(.data[[value_col]]),
      alpha2   = last(.data[[value_col]]),
      sex1     = first(Sex),
      sex2     = last(Sex),
      zygosity = first(zygosity),
      .groups  = "drop"
    )
  
  df_wide <- df_wide %>%
    mutate(
      sex1 = ifelse(sex1 %in% c("M","Male","1"), 1, 0),
      sex2 = ifelse(sex2 %in% c("M","Male","1"), 1, 0)
    )
  
  model <- umxACE(
    selDVs  = "alpha",
    selCovs = "sex",
    mzData  = df_wide %>% filter(zygosity %in% c("MZ","MZ/MZ/MZ")),
    dzData  = df_wide %>% filter(zygosity %in% c("DZ","OSDZ","DZ/OSDZ/DZ")),
    sep     = "",
    autoRun = TRUE 
  )
  
  return(list(data = df_wide, model = model))
}

#########################################################################
# Function 2: Extract Standardized Estimates, CIs, and P-values
#########################################################################
reportACE_clean <- function(model) {
  
  model <- mxModel(model, mxCI(c("top.A_std", "top.C_std", "top.E_std")))
  
  model_ci <- tryCatch(mxRun(model, intervals = TRUE), error = function(e) model)
  
  A_std <- mxEval(top.A_std[1,1], model_ci)
  C_std <- mxEval(top.C_std[1,1], model_ci)
  E_std <- mxEval(top.E_std[1,1], model_ci)
  
  ci_table <- summary(model_ci)$CI
  
  # Helper function to catch the CI regardless of how OpenMx formats the matrix name
  get_ci <- function(param_name, bound) {
    if(is.null(ci_table)) return(NA)
    
    base_name <- gsub("\\[1,1\\]", "", param_name)
    
    if(param_name %in% rownames(ci_table)) {
      return(ci_table[param_name, bound])
    } else if (base_name %in% rownames(ci_table)) {
      return(ci_table[base_name, bound])
    } else {
      return(NA)
    }
  }
  
  aeModel <- tryCatch(umxModify(model, update="c_r1c1", free=FALSE), error=function(e) NULL)
  ceModel <- tryCatch(umxModify(model, update="a_r1c1", free=FALSE), error=function(e) NULL)
  
  comp_p_A <- ifelse(!is.null(ceModel), umxCompare(model, ceModel)$p[2], NA)
  comp_p_C <- ifelse(!is.null(aeModel), umxCompare(model, aeModel)$p[2], NA)
  
  data.frame(
    Component = c("A", "C", "E"),
    Estimate  = round(c(A_std, C_std, E_std), 3),
    CI_lower  = round(c(
      get_ci("top.A_std[1,1]", "lbound"),
      get_ci("top.C_std[1,1]", "lbound"),
      get_ci("top.E_std[1,1]", "lbound")
    ), 3),
    CI_upper  = round(c(
      get_ci("top.A_std[1,1]", "ubound"),
      get_ci("top.C_std[1,1]", "ubound"),
      get_ci("top.E_std[1,1]", "ubound")
    ), 3),
    p_value   = c(comp_p_A, comp_p_C, NA) 
  )
}

#########################################################################
# Execution
#########################################################################
shannon <- read.csv("./inputs/shannon_comp_func_map.csv")

comp_T1 <- run_ace_timepoint(shannon, "T1", "shannon_composition") #shannon_pwy_stratified_filtered; # shannon_pwy_unstratified_withun
comp_T2 <- run_ace_timepoint(shannon, "T2", "shannon_composition")
comp_T3 <- run_ace_timepoint(shannon, "T3", "shannon_composition")

res_T1 <- reportACE_clean(comp_T1$model)
res_T2 <- reportACE_clean(comp_T2$model)
res_T3 <- reportACE_clean(comp_T3$model)




