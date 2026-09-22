##############################################################################
## data transformation; PCA plot
##############################################################################

library(dplyr)
library(tibble)
library(mixOmics)
library(factoextra)
library(ggplot2)

filter_low_abundance <- function(mat, min_pct = 0.01) {
  total_sum <- sum(mat, na.rm = TRUE)
  if (total_sum == 0) return(mat)
  keep_cols <- which((colSums(mat, na.rm = TRUE) * 100 / total_sum) > min_pct)
  return(mat[, keep_cols, drop = FALSE])
}

transform_clr <- function(count_df, filter_pct = 0.01) {
  mat <- as.matrix(count_df)
  
  if (!is.null(filter_pct) && filter_pct > 0) { 
    mat <- filter_low_abundance(mat, percent = filter_pct)
  }
  
  min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
  mat_offset  <- mat + (min_nonzero / 2)
  
  mat_tss <- mat_offset / rowSums(mat_offset)
  mat_clr <- mixOmics::logratio.transfo(mat_tss, logratio = "CLR")

  return(as.matrix(mat_clr))
}

clr_matrix <- transform_clr(raw_data, filter_pct = 0.01)

## check clustering on batch or timepoint
pca_res <- prcomp(clr_matrix, scale. = TRUE) #	

fviz_pca_ind(pca_res,label = "none",
  habillage = meta$Batch, addEllipses = TRUE, palette = c("red", "blue"),title = "bacterial composition", pointsize = 1) +
  theme_classic() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))

fviz_pca_ind(pca_res,label = "none", habillage =meta$Timepoint,addEllipses = TRUE, palette=c("T1"="#424242", "T2"="#E64626", "T3"="#0148A4"),title="bacterial composition", pointsize = 1) +
  theme_classic() + theme( panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.4) )

