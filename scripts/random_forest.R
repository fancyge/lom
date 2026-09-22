############################################################
## Random Forest classification: accuracy; error plot; VIP
############################################################
library(randomForest)
library(rsample)      
library(caret)        
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)

set.seed(123)
n_repeats <- 100
n_trees   <- 500
split_prop <- 0.75

# Preallocate collectors
err_list     <- vector("list", n_repeats)
vip_list     <- vector("list", n_repeats)
acc_records  <- numeric(n_repeats)

for (i in seq_len(n_repeats)) {
  
  data_split <- initial_split(data, prop = split_prop, strata = Timepoint)
  train <- training(data_split)
  test  <- testing(data_split)
  
  train$Timepoint <- factor(train$Timepoint)
  test$Timepoint  <- factor(test$Timepoint, levels = levels(train$Timepoint))
  
  rf <- randomForest(
    Timepoint ~ ., 
    data = train, 
    importance = TRUE, 
    ntree = n_trees
  )
  
  # tree convergence error rates
  err_df <- as.data.frame(rf$err.rate)
  err_df$Trees <- seq_len(nrow(err_df))
  err_df$Run   <- i
  err_list[[i]] <- err_df
  
  # Task B: Generalization performance (Test Set)
  pred <- predict(rf, test)
  cm   <- confusionMatrix(pred, test$Timepoint)
  acc_records[i] <- cm$overall["Accuracy"]
  
  # Task C: Variable Importance Measure (VIP)
  vip_raw <- as.data.frame(importance(rf)) %>%
    rownames_to_column(var = "Variable") %>%
    select(Variable, MeanDecreaseGini, MeanDecreaseAccuracy) %>%
    mutate(Run = i)
  vip_list[[i]] <- vip_raw
}

cat(sprintf(
  "Overall Accuracy across %d iterations: %.3f ± %.3f\n", 
  n_repeats, mean(acc_records), sd(acc_records)
))

####################################
### plot class-specific error rates
####################################
err_all <- bind_rows(err_list)

err_long <- err_all %>%
  pivot_longer(
    cols = -c(Trees, Run), 
    names_to = "Class", 
    values_to = "Error"
  ) %>%
  group_by(Trees, Class) %>%
  summarise(
    mean = mean(Error, na.rm = TRUE),
    sd   = sd(Error, na.rm = TRUE),
    .groups = "drop"
  )

class_colors <- c("OOB" = "dark green","T1"  = "#424242", "T2"  = "#E64626", "T3"  = "#0148A4")

p_err <- ggplot(err_long, aes(x = Trees, y = mean, color = Class, fill = Class)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.5, linetype = "twodash") +
  scale_color_manual(values = class_colors) +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "Random Forest Convergence & Out-of-Bag Error",
    x = "Number of Trees",
    y = "Error Rate"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title   = element_text(size = 10, face = "bold"),
    axis.line    = element_line(linewidth = 0.4, color = "black"),
    axis.text    = element_text(color = "black")
  )

print(p_err)

####################################
### variable importance
####################################
vip_combined <- bind_rows(vip_list)

vip_summary <- vip_combined %>%
  group_by(Variable) %>%
  summarise(
    mean_Gini     = mean(MeanDecreaseGini, na.rm = TRUE),
    sd_Gini       = sd(MeanDecreaseGini, na.rm = TRUE),
    mean_Accuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
    sd_Accuracy   = sd(MeanDecreaseAccuracy, na.rm = TRUE),
    n_models      = n(),
    .groups       = "drop"
  ) %>%
  arrange(desc(mean_Accuracy))

# ----------------------------------------------------------
####################################
##### Top 10 VIP & Abundance Distribution
####################################
top10_vip <- vip_summary %>% 
  slice_head(n = 10)

species_order <- top10_vip$Variable

# Dynamic expansion limits for lollipop plot
x_min <- floor(min(top10_vip$mean_Accuracy, na.rm = TRUE) * 0.9)
x_max <- ceiling(max(top10_vip$mean_Accuracy, na.rm = TRUE) * 1.15)

# Panel 1: Lollipop plot of Mean Decrease Accuracy or Mean Decrease Gini
p1 <- ggplot(top10_vip, aes(x = mean_Accuracy, y = factor(Variable, levels = rev(species_order)))) +
  geom_segment(
    aes(x = x_min, xend = mean_Accuracy, 
        y = factor(Variable, levels = rev(species_order)), 
        yend = factor(Variable, levels = rev(species_order))),
    color = "grey70", linewidth = 1.0
  ) +
  geom_point(size = 3, color = "#2B6CB0") +
  geom_text(aes(label = sprintf("%.2f", mean_Accuracy)), hjust = -0.3, size = 2.8) +
  scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Mean Decrease Accuracy", y = "Feature") +
  theme_minimal(base_size = 9) +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "black")
  )

# Panel 2: Abundance per Timepoint
mean_abund <- data %>%
  select(Timepoint, all_of(species_order)) %>%
  pivot_longer(-Timepoint, names_to = "Species", values_to = "Abundance") %>%
  group_by(Timepoint, Species) %>%
  summarise(MeanVal = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Species   = factor(Species, levels = rev(species_order)),
    Timepoint = factor(Timepoint, levels = c("T3", "T2", "T1"))
  )

p2 <- ggplot(mean_abund, aes(x = MeanVal, y = Species, fill = Timepoint)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("T1" = "#424242", "T2" = "#E64626", "T3" = "#0148A4")) +
  labs(x = "Mean Transformed Abundance", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(color = "black")
  )

# Combine panels side-by-side
figure_vip <- p1 + p2 + plot_layout(guides = "collect")
print(figure_vip)
