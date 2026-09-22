##############################################################################
#### ternary plot for Streptococcus
############################################################################## 

library(ggtern)
library(dplyr)

strep_ternary_data <- ace_species %>%
  filter(grepl("Streptococcus", feature)) %>%
  filter(!is.na(A) & !is.na(C) & !is.na(E))

p_dots <- ggtern(data = strep_ternary_data, 
                 aes(x = C, y = A, z = E)) + 
  
  geom_point(aes(color = Timepoint), size = 2.5, alpha = 0.9) +
  
  scale_color_manual(values = c("T1"="#424242", "T2"="#E64626", "T3"="#0148A4")) +
  
  labs(x = "C (Shared Env)", xarrow = "More C",
       y = "A (Genetics)", yarrow = "More A",
       z = "E (Unique Env)", zarrow = "More E") +
  
  theme_bw() + 
  theme_showarrows() + 
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    
    tern.axis.title.T = element_text(color = "#F9C013", face = "bold", size = 12),
    tern.axis.arrow.T = element_line(color = "#F9C013", linewidth = 1),
    tern.axis.text.T  = element_text(color = "#F9C013"),
    tern.axis.line.T = element_line(color = "#F9C013", linewidth = 0.5),
    tern.panel.grid.major.T = element_line(color = "#F9C013", linetype = "dashed", linewidth = 0.5),
    tern.axis.arrow.text.T = element_text(color = "#F9C013", face = "bold"),
    tern.axis.title.L = element_text(color = "#4BACC6", face = "bold", size = 12),
    tern.axis.arrow.L = element_line(color = "#4BACC6", linewidth = 1),
    tern.axis.text.L  = element_text(color = "#4BACC6"),
    tern.axis.line.L = element_line(color = "#4BACC6", linewidth = 0.5), 
    tern.panel.grid.major.L = element_line(color = "#4BACC6", linetype = "dashed", linewidth = 0.5),
    tern.axis.arrow.text.L = element_text(color = "#4BACC6", face = "bold"),

    tern.axis.title.R = element_text(color = "#1F4E78", face = "bold", size = 12),
    tern.axis.arrow.R = element_line(color = "#1F4E78", linewidth = 1),
    tern.axis.text.R  = element_text(color = "#1F4E78"),
    tern.axis.line.R = element_line(color = "#1F4E78", linewidth = 0.5), 
        tern.panel.grid.major.R = element_line(color = "#1F4E78", linetype = "dashed", linewidth = 0.5),
        tern.axis.arrow.text.R = element_text(color = "#1F4E78", face = "bold")
  )

p_dots


