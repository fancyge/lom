###################################################################
## MAG statistics for figures
###################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

theme_panel <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.text = element_text(color = "black", size = 10, face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.4)
  )


mags_df <- mags_all %>%
  mutate(
    Sample_ID = genome %>%
      str_remove("^(HQ_|MQ_)") %>%
      str_remove("_bin.*$"),
    Quality = case_when(
      completeness >= 90 & contamination < 5 ~ "HQ",
      completeness >= 50 & completeness < 90 & contamination < 5 ~ "MQ",
      TRUE ~ "Low"
    ),
    log_N50 = log10(N50 / 1000)
  ) %>%
  filter(Quality %in% c("HQ", "MQ"))

mags_meta <- mags_df %>%
  left_join(
    meta %>% mutate(Sample_ID = str_replace_all(Sequence_ID, "-", "_")) %>% select(Sample_ID, Timepoint),
    by = "Sample_ID"
  ) %>%
  filter(!is.na(Timepoint), Timepoint %in% c("T1", "T2", "T3"))

## 1. MAG Quality Distributions per Timepoint (T1, T2, T3)
long_quality_tp <- mags_meta %>%
  select(genome, Timepoint, completeness, contamination, log_N50) %>%
  pivot_longer(
    cols = c(completeness, contamination, log_N50),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      recode(Metric,
        completeness  = "% Completeness",
        contamination = "% Contamination",
        log_N50       = "log10(N50 in kbp)"
      ),
      levels = c("% Completeness", "% Contamination", "log10(N50 in kbp)")
    )
  )

p1_quality_tp <- ggplot(long_quality_tp, aes(x = Timepoint, y = Value, fill = Timepoint)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("T1" = "#424242", "T2" = "#E64626", "T3" = "#0148A4")) +
  labs(x = "Timepoint", y = NULL) +
  theme_panel +
  theme(legend.position = "none")

## 2. MAG quality stratified overall by quality (HQ vs MQ)
long_quality_tier <- mags_df %>%
  select(genome, Quality, completeness, contamination, log_N50) %>%
  pivot_longer(
    cols = c(completeness, contamination, log_N50),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      recode(Metric,
        completeness  = "% Completeness",
        contamination = "% Contamination",
        log_N50       = "log10(N50 in kbp)"
      ),
      levels = c("% Completeness", "% Contamination", "log10(N50 in kbp)")
    ),
    Quality = factor(Quality, levels = c("HQ", "MQ"))
  )

p2_quality_hq_mq <- ggplot(long_quality_tier, aes(x = Quality, y = Value, fill = Quality)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("HQ" = "#1b9e77", "MQ" = "#d95f02")) +
  labs(x = "Quality Tier", y = NULL) +
  theme_panel +
  theme(legend.position = "none")

## 3. Phylum composition pie chart
phylum_recode <- c(
  "Proteobacteria"             = "Pseudomonadota",
  "Firmicutes"                 = "Bacillota",
  "Actinobacteria"             = "Actinomycetota",
  "Bacteroidetes"              = "Bacteroidota",
  "Spirochaetes"               = "Spirochaetota",
  "Fusobacteria"               = "Fusobacteriota",
  "Mycoplasmatota"             = "others",
  "Tenericutes"                = "others",
  "Candidatus Gracilibacteria" = "others",
  "Thermodesulfobacteriota"    = "others"
)

phylum_colors <- c(
  "Pseudomonadota"             = "#E64626",
  "Bacillota"                  = "#0148A4",
  "Bacteroidota"               = "#FFB800",
  "Actinomycetota"             = "#007E3B",
  "Spirochaetota"              = "#4E98D3",
  "Candidatus Saccharibacteria" = "#FDCA90",
  "Campylobacterota"           = "#FBF38D",
  "Fusobacteriota"             = "#F79C72",
  "Deinococcota"               = "#BDDC96",
  "others"                     = "#757575"
)

phylum_df <- mags_phy %>%
  mutate(phylum_std = recode(phylum, !!!phylum_recode)) %>%
  count(phylum_std, name = "n") %>%
  mutate(
    pct = n / sum(n) * 100,
    label = ifelse(pct >= 3.0, sprintf("%.1f%%", pct), ""),
    phylum_std = reorder(phylum_std, -n)
  )

p3_phylum_pie <- ggplot(phylum_df, aes(x = "", y = n, fill = phylum_std)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
  scale_fill_manual(values = phylum_colors) +
  labs(fill = "Phylum") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

## 4. Stacked Bar Plot: Distribution of uSGBs by Phylum (FGBs vs GGBs)
usgb_df <- mags_phy %>%
  mutate(
    phylum_std = recode(phylum, !!!phylum_recode),
    Novelty_Tier = case_when(
      grepl("FGB|uFGB", GGB_FGB_status, ignore.case = TRUE) ~ "FGB (Family-level novelty)",
      grepl("GGB|uGGB", GGB_FGB_status, ignore.case = TRUE) ~ "GGB (Genus-level novelty)",
      TRUE ~ "GGB (Genus-level novelty)" # Fallback based on nomenclature pattern
    )
  ) %>%
  filter(grepl("uSGB|unknown|novel", SGB, ignore.case = TRUE)) %>%
  count(phylum_std, Novelty_Tier) %>%
  mutate(
    Novelty_Tier = factor(Novelty_Tier, levels = c("FGB (Family-level novelty)", "GGB (Genus-level novelty)"))
  )

p4_usgb_bar <- ggplot(usgb_df, aes(y = reorder(phylum_std, n, sum), x = n, fill = Novelty_Tier)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.65, color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "FGB (Family-level novelty)" = "#E64626", "GGB (Genus-level novelty)"  = "#0148A4"
    )
  ) +
  labs(
    y = "Phylum",
    x = "Number of uSGBs",
    fill = "Taxonomic Novelty"
  ) +
  theme_panel +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    axis.text.y = element_text(size = 9)
  )


print(p1_quality_tp)
print(p2_quality_hq_mq)
print(p3_phylum_pie)
print(p4_usgb_bar)
