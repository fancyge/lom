############################################################################################
#### run plsda to correct batch effect within function tables; using clr transformed data
############################################################################################
library(mixOmics)
library(PLSDAbatch)

func.trt.tune <- plsda(X = func, Y = grm_tp, ncomp = 5) 

func.trt.tune$prop_expl_var #choose the number that explains 100% variance in the outcome matrix Y,
sum(func.trt.tune$prop_expl_var$Y[seq_len(2)]) 


func.batch.tune <- PLSDA_batch(X = func, 
                             Y.trt = grm_tp, Y.bat = grm_batch,
                             ncomp.trt = 2, ncomp.bat = 10)

func.batch.tune$explained_variance.bat 
sum(func.batch.tune $explained_variance.bat$Y[seq_len(1)]) 

func.PLSDA_batch.res <- PLSDA_batch(X = func, 
                                  Y.trt = grm_tp, Y.bat = grm_batch,
                                  ncomp.trt = 2, ncomp.bat = 1)

# plsda batch corrected function table
func.PLSDA_batch.res$X.nobatch

## pca plot using corrected table
pca_res <- prcomp(func.PLSDA_batch.res$X.nobatch, scale. = TRUE) #	

fviz_pca_ind( pca_res,label = "none", habillage = grm_batch, addEllipses = TRUE, palette = c("red", "blue"),title = "pathway_unstratified", pointsize = 1) +
  theme_classic() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))

fviz_pca_ind(pca_res,label = "none", habillage = grm_tp,addEllipses = TRUE, palette=c("T1"="#424242", "T2"="#E64626", "T3"="#0148A4"),title="pathway_unstratified", pointsize = 1) + theme_classic() +  theme(  panel.grid = element_blank(), panel.border = element_rect(colour = "black",fill = NA,linewidth = 0.4))


############################################################################################
#### permannova to assess batch effect correction
############################################################################################

# before correction
adonis2(dist(as.matrix(clr_matrix) ~ grm_batch, method="euclidean") 
# after correction
adonis2(func.PLSDA_batch.res$X.nobatch ~ grm_batch, method="euclidean")











