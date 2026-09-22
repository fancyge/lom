### ordination plot using Nonmetric Multidimensional Scaling (NMDS) based on unifrac distance 

library(ggplot2)
library(vegan)

data <- read.delim("./inputs/unifrac_matrix.tsv")

set.seed(1)
nMDS <- metaMDS(data, distance = data, k = 2, maxit=30) 

scores <- as.data.frame(nMDS$points) 
meta <- read.delim("./inputs/metadata.txt", header = T, sep = "\t", row.names = 1) 
scores.meta <- merge(scores, meta, by = 'row.names')

scores.meta[,1] <- NULL

colnames(scores.meta) <- c("nMDS1", "nMDS2", names(meta)) 

dim(scores.meta)

ggplot(scores.meta, aes(nMDS1, nMDS2, color= Timepoint)) + 
  geom_point(size = 0.5, alpha = 1)  + 
  theme(axis.text = element_blank()) +ggtitle("Weighted UniFrac") + theme_bw() + stat_ellipse(show.legend = FALSE,linewidth = 0.2)+scale_colour_manual(values=c("T1"="#424242", "T2"="#E64626", "T3"="#0148A4"))+labs(color="Timepoints")+ theme(
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
       axis.line = element_blank(),  
    text = element_text(size = 10),
    legend.position = ""
  )

