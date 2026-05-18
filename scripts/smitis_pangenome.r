#############################################################################################
##### Streptococcus mitis ###################################################################
##### gene family richness; pangenome markers; UpSet plot; top20 variant gene families ######
#############################################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(vegan)
library(ggplot2)
library(broom)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(grid)

############################################################################################
######################## plot - gene family richness
############################################################################################
data1 <- read.delim("./inputs/smitis_panphlan_profile.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

meta_df <- read.delim("./inputs/metadata.txt", sep="\t", stringsAsFactors = FALSE) %>%
  select(sample = Sample_ID, timepoint = Timepoint) %>%
  filter(!is.na(timepoint) & timepoint %in% c("T1", "T2", "T3"))
  
family_ids <- data1[[1]]
sample_ids <- colnames(data1)[-1]

pa <- as.matrix(data1[ , -1])            
rownames(pa) <- family_ids
colnames(pa) <- sample_ids
pa_num <- (pa > 0) * 1

richness <- colSums(pa_num)               

diversity_df <- tibble(
  sample   = sample_ids,
  richness = richness
) %>% 
  left_join(meta_df, by = "sample") %>%
  filter(!is.na(timepoint))

set.seed(123)
min_n <- 60 #min(table(diversity_df$timepoint)) 
n_iter <- 1000

subsample_once <- function(df, min_n) {
  df %>%
    group_by(timepoint) %>%
    sample_n(min_n, replace = FALSE) %>%
    summarise(richness_mean = mean(richness), .groups = "drop")
}

subsample_results <- map_dfr(1:n_iter, ~{
  subsample_once(diversity_df, min_n) %>% mutate(iter = .x)
})

subsample_summary <- subsample_results %>%
  group_by(timepoint) %>%
  summarise(
    mean_richness = mean(richness_mean),
    sd_richness   = sd(richness_mean),
    ci_low        = quantile(richness_mean, 0.025),
    ci_high       = quantile(richness_mean, 0.975)
  )

p_gf_subsample <- ggplot(subsample_results, aes(x = timepoint, y = richness_mean, fill = timepoint)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, position = position_dodge(width = 0.2)) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 0.05, alpha = 0.7) +
  scale_fill_manual(name = "Timepoint", values = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")) +
  scale_color_manual(name = "Timepoint", values = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")) +
  theme_bw(base_size = 10) +coord_cartesian(ylim = c(1650, 1900))+
  theme(
    axis.title.x     = element_blank(),
    axis.text        = element_text(color = "black"),
    axis.ticks       = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "Gene family richness")

############################################################################################
######################## plot - number of pangenome markers
############################################################################################
map_stats <- read.delim("./inputs/smitis_panphlan_mapstats.txt", stringsAsFactors = FALSE)

map_stats2_filt <- map_stats %>%
  left_join(meta_df, by = "sample") %>%
  filter(!is.na(timepoint))

div_marker <- map_stats2_filt %>%
  select(sample, timepoint, richness = total_markers)

min_n2 <- 60 #min(table(div_marker$timepoint)) 

subsample_results_marker <- map_dfr(1:n_iter, ~{
  div_marker %>%
    group_by(timepoint) %>%
    sample_n(min_n2, replace = FALSE) %>%
    summarise(richness_mean = mean(richness), .groups = "drop") %>%
    mutate(iter = .x)
})

summary_marker <- subsample_results_marker %>%
  group_by(timepoint) %>%
  summarise(
    mean_richness = mean(richness_mean),
    sd_richness   = sd(richness_mean),
    ci_low        = quantile(richness_mean, 0.025),
    ci_high       = quantile(richness_mean, 0.975)
  )

p_marker_subsample <- ggplot(subsample_results_marker, aes(x = timepoint, y = richness_mean, fill = timepoint)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, position = position_dodge(width = 0.2), size = 0.5) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 0.05, alpha = 0.7) +
  scale_fill_manual(name = "Timepoint", values = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")) +
  scale_color_manual(name = "Timepoint", values = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")) +
  theme_bw(base_size = 10) +
  labs(
    y = "# Pangenome marker genes",
    x = NULL
  ) +coord_cartesian(ylim = c(50000, 90000))+
  theme(
    axis.title.y     = element_text(margin = ggplot2::margin(r = 10)),
    axis.text        = element_text(color = "black"),
    axis.ticks       = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
## combine the two plots
#library(patchwork)
combined_plot_subsample <- (p_gf_subsample + p_marker_subsample) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

############################################################################################
####################### UpSet plot - gene family distributions
############################################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(UpSetR)

pa_num <- (as.matrix(data1[, -1]) > 0) * 1
rownames(pa_num) <- data1[[1]]
sample_ids <- colnames(pa_num)

meta_df_filtered <- meta_df %>%
  filter(sample %in% sample_ids)

min_n <- 60 
n_iter <- 1000  

get_sets_iter <- function() {
  sampled_meta <- meta_df_filtered %>%
    group_by(timepoint) %>%
    sample_n(min_n, replace = FALSE)

  lapply(unique(sampled_meta$timepoint), function(tp) {
    samples_tp <- sampled_meta$sample[sampled_meta$timepoint == tp]
    fams_tp <- rownames(pa_num)[rowSums(pa_num[, samples_tp, drop = FALSE] > 0) > 0]
    fams_tp
  }) %>%
  setNames(unique(sampled_meta$timepoint))
}

set.seed(123) 
sets_list <- map(1:n_iter, ~get_sets_iter())

get_intersections <- function(sets, iter_index) {
  T1 <- sets[["T1"]]
  T2 <- sets[["T2"]]
  T3 <- sets[["T3"]]
  
  tibble(
    iter = iter_index,
    T1_only = length(setdiff(T1, union(T2, T3))),
    T2_only = length(setdiff(T2, union(T1, T3))),
    T3_only = length(setdiff(T3, union(T1, T2))),
    T1_T2 = length(setdiff(intersect(T1, T2), T3)),
    T1_T3 = length(setdiff(intersect(T1, T3), T2)),
    T2_T3 = length(setdiff(intersect(T2, T3), T1)),
    T1_T2_T3 = length(Reduce(intersect, list(T1, T2, T3)))
  )
}

intersections_df <- map_dfr(1:n_iter, function(i) {
  get_intersections(sets_list[[i]], i)
})

summary_intersections <- intersections_df %>%
  summarise(across(-iter, list(
    median = median,
    ci_low = ~quantile(.x, 0.025),
    ci_high = ~quantile(.x, 0.975)
  )))

exact_medians <- intersections_df %>% 
  summarise(across(-iter, median)) %>% 
  as.list()

# Format the input exactly how UpSetR expects it
expression_input <- c(
  "T1" = exact_medians$T1_only,
  "T2" = exact_medians$T2_only,
  "T3" = exact_medians$T3_only,
  "T1&T2" = exact_medians$T1_T2,
  "T1&T3" = exact_medians$T1_T3,
  "T2&T3" = exact_medians$T2_T3,
  "T1&T2&T3" = exact_medians$T1_T2_T3
)

# Plot the mathematically exact UpSet Plot
upset(
  fromExpression(expression_input),
  sets = c("T3", "T2", "T1"), 
  keep.order = TRUE,   
  main.bar.color = "#1F4E78",
  
  # Positional vector matching the sets order above (T3=Blue, T2=Orange, T1=Grey)
  sets.bar.color = c("#0148A4", "#E64626", "#424242"), 
  
  # Increased font sizes for publication readability
  text.scale = c(1.5, 1.2, 1.5, 1.2, 1.5, 1.5), 
  
  sets.x.label = "Median Set Size",
  mainbar.y.label = "Median Intersection Size",
  point.size = 2.5,
  line.size  = 1.0
)

c(C = "#4BACC6", A = "#F9C013", E = "#1F4E78")) 

############################################################################################
####################### pheatmap plot - top 20 gene family distributions
############################################################################################
keep <- rowSums(pa_num) > 5   # keep families present in >5 samples
pa_num_filt <- pa_num[keep, ]

meta_df <- meta_df[match(colnames(pa_num_filt), meta_df$sample), ]
stopifnot(all(colnames(pa_num_filt) == meta_df$sample))

test_gene <- function(presence, meta_df) {
  df <- data.frame(
    presence = presence,
    timepoint = factor(meta_df$timepoint)
  )
  fit <- try(glm(presence ~ timepoint, data = df, family = "binomial"), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    return(data.frame(term = NA, p.value = NA))
  }
  
  tidy(fit) %>%
    filter(term != "(Intercept)") %>%  # drop intercept
    summarise(p.value = min(p.value))  # take smallest p across contrasts
}

results <- apply(pa_num_filt, 1, function(x) test_gene(x, meta_df))
results_df <- bind_rows(results, .id = "gene_family")
results_df <- results_df %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p_adj)

sig_genes <- results_df %>% filter(!is.na(p_adj) & p_adj < 0.05)
top20 <- sig_genes$gene_family[1:20]

an <- read.delim("./inputs/top_gf_logistic.txt")
meta_df_filtered <- meta_df %>%
  filter(sample %in% sample_ids)

samples_keep <- meta_df_filtered %>%
  arrange(factor(timepoint, levels = c("T1", "T2", "T3"))) %>%
  pull(sample)

mat <- pa_num[top20, samples_keep]

col_ha <- HeatmapAnnotation(
  Timepoint = factor(meta_df_filtered $timepoint[match(colnames(mat), meta_df_filtered $sample)],
                     levels = c("T1", "T2", "T3")),
  col = list(Timepoint = c(T1 = "#424242", T2 = "#E64626", T3 = "#0148A4")),
   annotation_name_gp = gpar(fontsize = 0),
  annotation_name_side = "left"
)

matched_annotations <- an$annotation[match(rownames(mat), an$NR90)]
labels_for_plot <- paste(gsub("UniRef90_","",rownames(mat)), matched_annotations)

row_ha <- rowAnnotation(
  Annotation = anno_text(
    labels_for_plot,
    gp = gpar(fontsize = 9),
    just = "left"
  )
)

hc_col <- hclust(dist(t(mat)), method = "complete")
dend_col <- as.dendrogram(hc_col)

dend_col[[2]] <- rev(dend_col[[2]])
hc_col2 <- as.hclust(dend_col)

ht <- Heatmap(
  mat,
  name = "Presence",
  col = c("0" = "#A7CEE1", "1" = "#00A485"),
  cluster_rows = TRUE,
  cluster_columns = hc_col, 
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_rot = 0,  
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
