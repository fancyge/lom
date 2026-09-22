#############################################################################################
########## Strain-level analysis on Streptococcus mitis 
########## Evaluating the influence of timepoint on pairwise genetic distance 
########## via PERMANOVA statistics and PCoA ordination
#############################################################################################
library(dplyr) 
library(vegan)
library(stringr)
library(ggplot2)

data <- read.table("./inputs/smitis_distmat.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.delim("./inputs/metadata.txt", sep="\t", stringsAsFactors = FALSE)

meta_aligned <- data.frame(SampleID = rownames(data)) %>%
  left_join(metadata, by = c("SampleID" = "Sample_ID")) %>%
  select(SampleID, Family_ID, tp = Timepoint, Sex=Sex) %>%
  mutate(
    tp = ifelse(grepl("GCF_", SampleID), "REF", tp),
    Family_ID = factor(Family_ID),
    Sex=factor(Sex)
  )

### statistical tests
biological_samples <- meta_aligned %>% filter(tp != "REF" & !is.na(tp)) %>% pull(SampleID)

data_stats <- data[biological_samples, biological_samples]
meta_stats <- meta_aligned %>% filter(SampleID %in% biological_samples) %>% mutate(Family_ID = factor(Family_ID),Sex=factor(Sex))

dist_mat_stats <- as.dist(data_stats)
adonis_result <- adonis2(dist_mat_stats ~ tp+Sex, data = meta_stats, strata = meta_stats$Family_ID, permutations = 1000)
print(adonis_result)

bd <- betadisper(dist_mat_stats, meta_stats$tp)
anova(bd)
permutest(bd)

### PCoA plotting

dist_mat_full <- as.dist(data)
e.sir.pcoa <- cmdscale(dist_mat_full, eig = TRUE)

variance <- head(eigenvals(e.sir.pcoa) / sum(eigenvals(e.sir.pcoa)))
x_variance <- round(variance[1] * 100, 1)
y_variance <- round(variance[2] * 100, 1)

scores <- as.data.frame(e.sir.pcoa$points)
colnames(scores) <- c("PCo1", "PCo2")
scores$SampleID <- rownames(scores)

scores_meta <- left_join(scores, meta_aligned, by = "SampleID")

pcoa_plot <- ggplot(scores_meta, aes(x = PCo1, y = PCo2, color = tp)) + 
  geom_point(size = 1.5, alpha = 0.75)  + 
  # Only draw ellipses for the actual biological groups, not the individual REF points
  stat_ellipse(data = subset(scores_meta, tp != "REF"), 
               aes(group = tp), level = 0.95, linetype = "dashed") +
  scale_color_manual(
    name   = "Timepoint",
    values = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4", REF = "#007E3B")
  ) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 8, color = "black"),
    axis.title        = element_text(size = 10),
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    axis.ticks.length = unit(1, "pt"),       
    axis.ticks        = element_line(linewidth = 0.2) 
  ) +
  xlab(paste0("PCo1 (", x_variance, "% variance explained)")) +
  ylab(paste0("PCo2 (", y_variance, "% variance explained)")) 

print(pcoa_plot)

