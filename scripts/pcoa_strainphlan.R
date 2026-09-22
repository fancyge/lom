################################################################################
## ordination plot using the distance matrix calculated from StrainPhlAn msa
################################################################################

library(ggplot2)
library(vegan)

plot_strainphlan_pcoa <- function(species, base_dir, output_file) {

  dir <- file.path(base_dir, paste0("output_", species))

  data <- read.table(file.path(dir, "distmat.txt"), skip = 8, fill = TRUE,
                     sep = "\t", strip.white = TRUE, check.names = FALSE)

  samples <- gsub(" .*", "", data[[ncol(data)]])

  data[1] <- NULL
  data[ncol(data)] <- NULL
  data[ncol(data)] <- NULL

  rownames(data) <- samples
  colnames(data) <- samples
  data[lower.tri(data)] <- t(data)[lower.tri(data)]

  pcoa <- cmdscale(data, eig = TRUE)
  variance <- head(eigenvals(pcoa) / sum(eigenvals(pcoa)))

  scores <- as.data.frame(pcoa$points)
  meta <- read.delim(file.path(dir, "tp.txt"), header = TRUE,
                     sep = "\t", row.names = 1)

  scores <- merge(scores, meta, by = "row.names")
  rownames(scores) <- scores[, 1]
  scores[, 1] <- NULL
  colnames(scores) <- c("PCo1", "PCo2", "sample_type")

  ggplot(scores, aes(PCo1, PCo2, color = sample_type)) +
    geom_point(size = 1.2, alpha = 0.75) +
    stat_ellipse(aes(group = sample_type), level = 0.95, linetype = "dashed") +
    scale_color_manual(
      name = "Timepoint",
      values = c(T1 = "#424242", T2 = "#E64626",
                 T3 = "#0148A4", REF = "#007E3B")
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.ticks.length = unit(1, "pt"),
      axis.ticks = element_line(size = 0.2),
      plot.title = element_text(size = 12)
    ) +
    labs(
      x = paste0("PCo1 (", as.integer(variance[1] * 100), "% variance explained)"),
      y = paste0("PCo2 (", as.integer(variance[2] * 100), "% variance explained)"),
      title = species
    ) -> p

  ggsave(output_file, p, width = 6, height = 5)
}

