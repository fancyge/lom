###################################################################################
##### alpha diversity of timepoint groups
##### shannon plot; statistics
###################################################################################
library(vegan)
library(ggplot2)
library(vegan)
library(lmerTest)
library(emmeans)


df <- read.csv("./inputs/shannon_comp_func_map.csv")
df$Timepoint <- factor(df$Timepoint, levels = c("T1", "T2", "T3"))

age_positions <- df %>%
  group_by(Timepoint) %>%
  summarise(mean_age = mean(Age.months., na.rm = TRUE))

df <- df %>%
  left_join(age_positions, by = "Timepoint")

############################## Bacterial composition ######################################## 
### plots 
ggplot(df, aes(x = Age.months., y = shannon_composition)) +
  geom_violin(aes(x = mean_age, group = Timepoint, fill = Timepoint),
              width = 20, alpha = 0.9, color = NA, trim = FALSE) +
  geom_point(aes(color = Timepoint), alpha = 0.75,size=0.3) +
  geom_smooth(color = "black", se = TRUE,size=0.7) +guides(
    color = guide_legend(override.aes = list(shape = 16, size = 2, fill = NA)),
    fill = guide_legend(override.aes = list(shape = 16, size = 2))
  )+
  scale_color_manual(values = c("#424242", "#E64626", "#0148A4")) +
  scale_fill_manual(values = c("#424242", "#E64626", "#0148A4")) +
  labs(
  title="Bacterial Composition",
    y = "Shannon Index",
    x = "Age (months)",
    fill = "",
    color = ""
  ) +theme_bw() +
  theme(
    panel.border = element_blank(),  # Remove full border
    text = element_text(size = 10),
    axis.line = element_line(color = "black"), 
    axis.line = element_line(color = "black",size=0.3),  
    legend.position = ""
  )

### stats 
lmer.shannon.comp.Time <- lmer(shannon_composition ~ Timepoint + Sex + (1|Family_ID) , data=df)
summary(lmer.shannon.comp.Time)

aov.shannon.comp.time <- anova(lmer.shannon.comp.Time)
summary(aov.shannon.comp.time)
aov.shannon.comp.time$Pr

contrasts_Shannon_comp_time <- emmeans(lmer.shannon.comp.Time, pairwise ~ Timepoint)
contrasts_Shannon_comp_time$contrasts


############################## Functional Diversity using filtered stratified pathways ################################# 
### plots 
ggplot(df, aes(x = Age.months., y = shannon_pwy_stratified_filtered)) +
  geom_violin(aes(x = mean_age, group = Timepoint, fill = Timepoint),
              width = 20, alpha = 0.9, color = NA, trim = FALSE) +
  geom_point(aes(color = Timepoint), alpha = 0.75,size=0.3) +
  geom_smooth(color = "black", se = TRUE,size=0.7) +guides(
    color = guide_legend(override.aes = list(shape = 16, size = 2, fill = NA)),
    fill = guide_legend(override.aes = list(shape = 16, size = 2))
  )+
  scale_color_manual(values = c("#424242", "#E64626", "#0148A4")) +
  scale_fill_manual(values = c("#424242", "#E64626", "#0148A4")) +
  labs(
  title="Functional Diversity",
    y = "Shannon Index",
    x = "Age (months)",
    fill = "",
    color = ""
  ) +theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),  # Remove full border
    text = element_text(size = 10),
    axis.line = element_line(color = "black"),  
    axis.line = element_line(color = "black",size=0.3),
    legend.position = ""
  )

### stats 
lmer.shannon.func.Time <- lmer(shannon_pwy_stratified_filtered ~ Timepoint + Sex + (1|Family_ID) , data=df)
summary(lmer.shannon.func.Time)

aov.shannon.func.time <- anova(lmer.shannon.func.Time)
summary(aov.shannon.func.time)
aov.shannon.func.time$Pr

contrasts_Shannon_func_time <- emmeans(lmer.shannon.func.Time, pairwise ~ Timepoint)
contrasts_Shannon_func_time$contrasts

