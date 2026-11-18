################################################################################
## Pairwise phylogenetic distances of S. mitis across developmental stages and relationship types
## Node-label permutation tests accounting for overlapping pairwise distances
################################################################################

library(ape)
library(dplyr)
library(ggplot2)
library(ggpubr)

set.seed(123)
n_perm <- 9999

tree <- read.tree("RAxML_bestTree.s__Streptococcus_mitis.StrainPhlAn3.tre")

dist_matrix <- cophenetic.phylo(tree)

meta_use <- meta %>%
  transmute(Sample_ID = as.character(Sample_ID),
            Twin_ID = as.character(Twin_ID),
            Family_ID = as.character(Family_ID),
            Timepoint = as.character(Timepoint)) %>%
  filter(Timepoint %in% c("T1", "T2", "T3"),
         !is.na(Twin_ID), !is.na(Family_ID)) %>%
  distinct()

sample_ids <- intersect(rownames(dist_matrix), meta_use$Sample_ID)
meta_use <- meta_use[match(sample_ids, meta_use$Sample_ID), ]
dist_matrix <- dist_matrix[sample_ids, sample_ids]

pair_idx <- which(lower.tri(dist_matrix), arr.ind = TRUE)

i1 <- pair_idx[, 2]
i2 <- pair_idx[, 1]

# Extract pairwise distances and annotate metadata
dist_annotated <- data.frame(
  Sample1 = sample_ids[i1],
  Sample2 = sample_ids[i2],
  Distance = dist_matrix[cbind(i1, i2)],
  Twin_ID1 = meta_use$Twin_ID[i1],
  Twin_ID2 = meta_use$Twin_ID[i2],
  Family_ID1 = meta_use$Family_ID[i1],
  Family_ID2 = meta_use$Family_ID[i2],
  Time1 = meta_use$Timepoint[i1],
  Time2 = meta_use$Timepoint[i2]
) %>%
  mutate(
    Time_order1 = match(Time1, c("T1", "T2", "T3")),
    Time_order2 = match(Time2, c("T1", "T2", "T3")),
    Time_low = c("T1", "T2", "T3")[pmin(Time_order1, Time_order2)],
    Time_high = c("T1", "T2", "T3")[pmax(Time_order1, Time_order2)],
    Time_pair = paste(Time_low, Time_high, sep = "-"),

    Comparison_Group = case_when(
      Twin_ID1 == Twin_ID2 & Time1 != Time2 & Time_pair == "T1-T2" ~ "Longitudinal (T1-T2)",

      Twin_ID1 == Twin_ID2 & Time1 != Time2 & Time_pair == "T2-T3" ~ "Longitudinal (T2-T3)",

      Twin_ID1 == Twin_ID2 & Time1 != Time2 & Time_pair == "T1-T3" ~ "Longitudinal (T1-T3)",

      Twin_ID1 != Twin_ID2 & Family_ID1 == Family_ID2 & Time1 == Time2 & Time1 == "T1" ~ "Within-family (T1)",

      Twin_ID1 != Twin_ID2 & Family_ID1 == Family_ID2 & Time1 == Time2 & Time1 == "T2" ~ "Within-family (T2)",

      Twin_ID1 != Twin_ID2 & Family_ID1 == Family_ID2 & Time1 == Time2 & Time1 == "T3" ~ "Within-family (T3)",

      Family_ID1 != Family_ID2 ~ "Unrelated",

      TRUE ~ "Other_Related"
    )
  ) %>%
  filter(Comparison_Group != "Other_Related")

# Check comparison numbers
print(table(dist_annotated$Comparison_Group))

# Define six comparisons against timepoint-matched unrelated pairs
test_groups <- c(
  "Within-family (T1)",
  "Within-family (T2)",
  "Within-family (T3)",
  "Longitudinal (T1-T2)",
  "Longitudinal (T2-T3)",
  "Longitudinal (T1-T3)"
)

# Define the matching developmental-stage combination for each test
test_time_pairs <- c("T1-T1","T2-T2","T3-T3","T1-T2","T2-T3","T1-T3")

# Observed Wilcoxon rank-sum statistics
rank_sum_stat <- function(x, y) {
  ranks <- rank(c(x, y), ties.method = "average")
  sum(ranks[seq_along(x)])
}

obs_W <- numeric(length(test_groups))
names(obs_W) <- test_groups

test_indices <- vector("list", length(test_groups))
names(test_indices) <- test_groups

for (g in test_groups) {

  tp_pair <- test_time_pairs[match(g, test_groups)]

  x_idx <- which(dist_annotated$Comparison_Group == g)

  y_idx <- which(
    dist_annotated$Comparison_Group == "Unrelated" &
      dist_annotated$Time_pair == tp_pair
  )

  if (length(x_idx) == 0 || length(y_idx) == 0) {
    stop(paste("Insufficient observations for", g))
  }

  test_indices[[g]] <- list(x = x_idx, y = y_idx)

  obs_W[g] <- rank_sum_stat(
    dist_annotated$Distance[x_idx],
    dist_annotated$Distance[y_idx]
  )
}

# Permutation test
# Shuffle sample labels within each developmental stage. The distance matrix remains fixed, preserving its pairwise network structure.
perm_W <- matrix(
  NA_real_,
  nrow = n_perm,
  ncol = length(test_groups),
  dimnames = list(NULL, test_groups)
)

n <- length(sample_ids)

stage_indices <- split(
  seq_len(n),
  factor(meta_use$Timepoint, levels = c("T1", "T2", "T3"))
)

for (b in seq_len(n_perm)) {

  if (b %% 20 == 0) message("Permutation ", b, " / ", n_perm)
  # Permute sample labels independently within each timepoint
  perm_map <- seq_len(n)
  
  for (tp in names(stage_indices)) {
    idx <- stage_indices[[tp]]
    perm_map[idx] <- sample(idx, length(idx))
  }
  # Recalculate permuted pairwise distances using the original matrix
  perm_dist <- dist_matrix[
    cbind(perm_map[i1], perm_map[i2])
  ]
  # Calculate all six permuted statistics
  for (g in test_groups) {
    x_idx <- test_indices[[g]]$x
    y_idx <- test_indices[[g]]$y
    perm_W[b, g] <- rank_sum_stat(
      perm_dist[x_idx],
      perm_dist[y_idx]
    )
  }
}

# Calculate two-sided empirical P values
empirical_results <- data.frame(
  Comparison = test_groups,
  group1 = "Unrelated",
  group2 = test_groups,
  Observed_W = as.numeric(obs_W),
  Empirical_P = NA_real_
)

for (g in test_groups) {
  null_mean <- mean(perm_W[, g])
  empirical_results$Empirical_P[
    empirical_results$Comparison == g
  ] <- (
    1 + sum(abs(perm_W[, g] - null_mean) >=
              abs(obs_W[g] - null_mean))
  ) / (n_perm + 1)
}

# Holm correction across six comparisons
empirical_results <- empirical_results %>%
  mutate(
    Adjusted_P = p.adjust(Empirical_P, method = "holm"),
    p.signif = case_when(
      Adjusted_P <= 0.0001 ~ "****",
      Adjusted_P <= 0.001 ~ "***",
      Adjusted_P <= 0.01 ~ "**",
      Adjusted_P <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(empirical_results)

plot_levels <- c(
  "Unrelated",
  "Within-family (T1)",
  "Within-family (T2)",
  "Within-family (T3)",
  "Longitudinal (T1-T2)",
  "Longitudinal (T2-T3)",
  "Longitudinal (T1-T3)"
)

plot_data <- dist_annotated %>%
  filter(Comparison_Group %in% plot_levels) %>%
  mutate(
    Comparison_Group = factor(
      Comparison_Group,
      levels = plot_levels
    )
  )

max_y <- max(plot_data$Distance, na.rm = TRUE)
step <- max_y * 0.08

stat_test <- empirical_results %>%
  mutate(
    y.position = max_y + step * row_number(),
    p_label = p.signif
  )

# plot
p <- ggplot(
  plot_data,
  aes(x = Comparison_Group, y = Distance, fill = Comparison_Group)
) +
  geom_boxplot(
    outlier.alpha = 0.1,
    outlier.size = 0.5,
    width = 0.6
  ) +
  stat_pvalue_manual(
    stat_test,
    label = "p_label",
    tip.length = 0.01,
    step.increase = 0,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  theme_minimal(base_size = 12) +
  labs(
    y = "Pairwise Phylogenetic Distance (Cophenetic)",
    x = ""
  ) +
  theme(
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5
    ),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text = element_text(color = "black", size = 9)
  )

print(p)
