################################################################################
## replication rates visualizations; group and individual species plots across Timepoints
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(forcats)

data <- read.delim("./ReplicationRates_combined_7groups.txt", stringsAsFactors = FALSE)

group_colors <- c("Actinomyces" = "#424242", "Corynebacterium" = "#E64626", 
                  "Fusobacterium" = "#0148A4", "Haemophilus" = "#FFB800", 
                  "Neisseria" = "#007E3B", "Streptococcus" = "#4E98D3", "TM7" = "#F79C72")

tp_colors <- c("T1" = "#424242", "T2" = "#E64626", "T3" = "#0148A4")

data <- data %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3")),
    group = factor(group, levels = rev(c(
      "Fusobacterium", "Streptococcus", "Corynebacterium",
      "Haemophilus", "Neisseria", "Actinomyces", "TM7"
    )))
  )

theme_irep <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = 11),
    text             = element_text(color = "black"),
    axis.text        = element_text(color = "black", size = 9),
    strip.text       = element_text(color = "black", size = 10),
    legend.text      = element_text(color = "black"),
    legend.title     = element_text(color = "black")
  )


## density plot
plot_ridges <- ggplot(data, aes(x = index, y = group, fill = Timepoint)) +
  geom_density_ridges(alpha = 0.6, scale = 1.3) +
  scale_fill_manual(values = tp_colors) +
  labs(title = "Density of Replication Indices", x = "Replication Index", y = "Group") +
  theme_irep + theme(legend.position = "right")

print(plot_ridges)


## plot for individual species
plot_single_species <- function(target_species) {
  
  test_data <- data %>% filter(species == target_species)
  
  if(nrow(test_data) == 0) {
    stop("Species not found in the dataset. Please check your spelling.")
  }
  
  p <- ggplot(test_data, aes(x = Timepoint, y = index)) +
    geom_boxplot(fill = "skyblue", width = 0.5, alpha = 0.5, outlier.shape = NA) +  
    geom_jitter(width = 0.12, alpha = 0.7, size = 1.5, color = "darkorange") +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 2) +
    labs(title = target_species, x = "Timepoint", y = "Replication index") +
    theme_irep + 
    theme(plot.title = element_text(face = "italic"))
  
  return(p)
}

plot_single_species("Fuso_pseudoperiodonticum")
plot_single_species("Neis_mucosa")
