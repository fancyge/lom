##############################################################################
#### top ACE features per tp 
############################################################################## 

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext) 

plot_unified_panel_c <- function(data, n_top = 10, y_label_face = "italic") {
  
  top_C_T1 <- data %>% filter(Timepoint == "T1", !is.na(C)) %>% arrange(desc(C)) %>% slice_head(n = n_top) %>% mutate(Group = "High Shared Env (C) - T1", SortVal = C)
  top_C_T2 <- data %>% filter(Timepoint == "T2", !is.na(C)) %>% arrange(desc(C)) %>% slice_head(n = n_top) %>% mutate(Group = "High Shared Env (C) - T2", SortVal = C)
  top_A_T3 <- data %>% filter(Timepoint == "T3", !is.na(A)) %>% arrange(desc(A)) %>% slice_head(n = n_top) %>% mutate(Group = "High Genetics (A) - T3", SortVal = A)
  
  panel_data <- bind_rows(top_C_T1, top_C_T2, top_A_T3) %>%
    pivot_longer(cols = c(A, C, E), names_to = "component", values_to = "variance") %>%
    mutate(
      Group = factor(Group, levels = c("High Shared Env (C) - T1", "High Shared Env (C) - T2", "High Genetics (A) - T3")),
      feature = tidytext::reorder_within(feature, SortVal, Group)
    )
  
  panel_data$component <- factor(panel_data$component, levels = c("E", "C", "A"))
  
  p <- ggplot(panel_data, aes(x = feature, y = variance, fill = component)) +
    geom_bar(stat = "identity", width = 0.5) + 
    scale_fill_manual(values = c(C = "#4BACC6", A = "#F9C013", E = "#1F4E78")) +
    coord_flip() +
    
    tidytext::scale_x_reordered() + 
    
    facet_wrap(~ Group, scales = "free_y", nrow = 3) + 
    
    labs(x = NULL, y = "Variance explained (%)", fill = "Component") +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 11,  color = "black", margin = margin(b=10)),
      axis.text.y = element_text(face = y_label_face, color = "black", size = 8),
      axis.text.x = element_text(color = "black", size = 8),
      axis.title.x = element_text(color = "black", size = 10, margin = margin(t=10)),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm")
    )
  
  return(p)
}

panel_C_species <- plot_unified_panel_c(data = ace_species, n_top = 10, y_label_face = "italic")
panel_C_pathways <- plot_unified_panel_c(data = ace_pwys, n_top = 10, y_label_face = "plain")



