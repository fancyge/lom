##############################################################################
### ACE bar plot 
############################################################################## 
library(ggplot2)
library(tidyr)
library(dplyr)

 ace <- read.csv("./inputs/ACE.csv")

ace_long <- ace %>%
  pivot_longer(
    cols = c("A", "C", "E"),
    names_to = "Component",
    values_to = "Value" ) %>%
  
  mutate(
    Time = sub(".*_", "", cat),
    Type = sub("_[^_]+$", "", cat)) %>% 
  
  group_by(cat) %>%
  mutate(  Percent = Value / sum(Value) * 100) %>%
  
  ungroup() %>% mutate(
    Type = factor( Type,  levels = c(  "comp", "pwy_strat",  "pwy_unstrat"),
    labels = c("Species","Species-stratified pathways", "Unstratified pathways" ) ),
    Time = factor(Time,levels = c("T1", "T2", "T3")  ), 
    Component = factor( Component, levels = c("A", "C", "E") ) )
  
  
ggplot(
  ace_long,
  aes(
    x = Time,
    y = Percent,
    fill = Component )) +  geom_col(
    position = position_dodge(width = 0.65), width = 0.75, colour = "white" ) +  
    geom_text(
    aes(label = round(Percent, 1)),
    position = position_dodge(width = 0.75), vjust = -0.3, size = 2.8  ) +  facet_wrap(  ~ Type, nrow = 1  ) +  scale_fill_manual(
    values = c( C = "#4BACC6", A = "#F9C013",E = "#1F4E78"    )) +  scale_y_continuous( limits = c(0, 100),expand = expansion(mult = c(0, 0.08))) +  labs(
    x = NULL,
    y = "Proportion of variance (%)",
    fill = "Variance component") +  theme_classic() +  theme(
    strip.background = element_blank(),    
    # Frame around each panel
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5 ),    
    legend.position = "bottom" )
  