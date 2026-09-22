
#########################################################################################
## panphlan plot for s. mitis: pangenome gene markers richness; upset plot; 
#########################################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(UpSetR)

set.seed(123)
n_iter <- 1000

### gene markers richness
meta_df <- read.delim("./Panphlan_Streptococcus_mitis/profile_results_920/meta_296.txt",stringsAsFactors = FALSE) %>%
  rename(timepoint = Timepoint) %>%
  mutate(sample = sample %>% str_remove("^result_map_") %>% str_remove("_smitis\\.csv$")) %>%
  filter(timepoint %in% c("T1", "T2", "T3")) %>%
  select(sample, timepoint)

data1 <- read.delim(
  "./Panphlan_Streptococcus_mitis/profile_results_920/result_profile_smitis_no_ref.tsv",
  stringsAsFactors = FALSE
)

pa_mat <- as.matrix(data1[, -1]) > 0
rownames(pa_mat) <- data1[[1]]
storage.mode(pa_mat) <- "numeric"

df_gf <- tibble(
  sample = colnames(pa_mat),
  richness = colSums(pa_mat)
) %>%
  inner_join(meta_df, by = "sample")

min_n_gf <- 60

subsample_results_gf <- map_dfr(seq_len(n_iter), function(i) {
  df_gf %>%
    group_by(timepoint) %>%
    slice_sample(n = min_n_gf, replace = FALSE) %>%
    summarise(richness_mean = mean(richness), .groups = "drop") %>%
    mutate(iter = i)
})

map_stats <- read.delim(
  "./panphlan-map-stats.txt",
  stringsAsFactors = FALSE
)

df_marker <- map_stats %>%
  inner_join(meta_df, by = "sample") %>%
  select(sample, timepoint, richness = total_markers)

min_n_marker <- 60

subsample_results_marker <- map_dfr(seq_len(n_iter), function(i) {
  df_marker %>%
    group_by(timepoint) %>%
    slice_sample(n = min_n_marker, replace = FALSE) %>%
    summarise(richness_mean = mean(richness), .groups = "drop") %>%
    mutate(iter = i)
})

tp_colors <- c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")

theme_clean <- theme_bw(base_size = 9) +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_gf_subsample <- ggplot(subsample_results_gf, aes(x = timepoint, y = richness_mean, fill = timepoint)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 0.05, alpha = 0.7) +
  scale_fill_manual(name = "Timepoint", values = tp_colors) +
  scale_color_manual(name = "Timepoint", values = tp_colors) +
  labs(y = "Gene family richness") +
  theme_clean

p_marker_subsample <- ggplot(subsample_results_marker, aes(x = timepoint, y = richness_mean, fill = timepoint)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 0.05, alpha = 0.7) +
  scale_fill_manual(name = "Timepoint", values = tp_colors) +
  scale_color_manual(name = "Timepoint", values = tp_colors) +
  labs(y = "# Pangenome marker genes") +
  theme_clean

combined_plot <- (p_gf_subsample + p_marker_subsample) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)


### upset plot
meta_upset <- meta_df %>%
  filter(sample %in% colnames(pa_mat))

min_n_upset <- 60

sets_list <- replicate(n_iter, {
  sampled <- meta_upset %>%
    group_by(timepoint) %>%
    slice_sample(n = min_n_upset, replace = FALSE) %>%
    ungroup()
  
  split_samples <- split(sampled$sample, sampled$timepoint)
  lapply(split_samples, function(s) {
    rownames(pa_mat)[rowSums(pa_mat[, s, drop = FALSE]) > 0]
  })
}, simplify = FALSE)

intersections_df <- map_dfr(seq_len(n_iter), function(i) {
  sets <- sets_list[[i]]
  t1 <- sets[["T1"]]
  t2 <- sets[["T2"]]
  t3 <- sets[["T3"]]
  
  tibble(
    iter = i,
    T1_only = length(setdiff(t1, union(t2, t3))),
    T2_only = length(setdiff(t2, union(t1, t3))),
    T3_only = length(setdiff(t3, union(t1, t2))),
    T1_T2 = length(setdiff(intersect(t1, t2), t3)),
    T1_T3 = length(setdiff(intersect(t1, t3), t2)),
    T2_T3 = length(setdiff(intersect(t2, t3), t1)),
    T1_T2_T3 = length(intersect(intersect(t1, t2), t3))
  )
})

median_iter <- intersections_df %>%
  mutate(diff = abs(T1_T2_T3 - median(T1_T2_T3))) %>%
  arrange(diff) %>%
  slice(1) %>%
  pull(iter)

upset(
  fromList(sets_list[[median_iter]]),
  sets = c("T3", "T2", "T1"),
  keep.order = TRUE,
  main.bar.color = "#007E3B",
  sets.bar.color = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4"),
  text.scale = c(1.2, 0.8, 1.0, 0.8, 1.2, 0.8),
  sets.x.label = "Set Size",
  point.size = 2.3,
  line.size = 0.7
)


### heatmap for top gf with annotation
an <- read.delim("./Panphlan_Streptococcus_mitis/top_gf_logistic.txt")

samples_keep <- meta_df %>%
  arrange(factor(timepoint, levels = c("T1", "T2", "T3"))) %>%
  pull(sample)

top20 <- an$NR90[1:20]

mat <- pa_num[top20, samples_keep]

col_ha <- HeatmapAnnotation(
  Timepoint = factor(meta_df$timepoint[match(colnames(mat), meta_df$sample)],
                     levels = c("T1", "T2", "T3")),
  col = list(Timepoint = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")),
   annotation_name_gp = gpar(fontsize = 0),
  annotation_name_side = "left"
)

top20_annot <- an %>%
  filter(NR90 %in% top20) %>%
  mutate(label = paste(NR90, annotation)) %>%
  arrange(match(NR90, rownames(mat)))  # ensure same order

row_ha <- rowAnnotation(
  Annotation = anno_text(
    top20_annot$label,
    gp = gpar(fontsize = 8),
    just = "left"
  )
)

hc_col <- hclust(dist(t(mat)), method = "complete")
dend_col <- as.dendrogram(hc_col)

dend_col[[2]] <- rev(dend_col[[2]])
hc_col2 <- as.hclust(dend_col)

a <- usyd_palette("extended")

ht <- Heatmap(
  mat,
  name = "Presence",
  col = c("0" = "#A7CEE1", "1" = a[11]),
  cluster_rows = TRUE,
  cluster_columns = hc_col, #hc_col2
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_rot = 0,  # horizontal labels
  column_names_gp = gpar(fontsize = 3),
  top_annotation = col_ha,
  right_annotation = row_ha,
  row_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("absent", "present"),
    title_position = "topcenter"
  )
)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE)




















