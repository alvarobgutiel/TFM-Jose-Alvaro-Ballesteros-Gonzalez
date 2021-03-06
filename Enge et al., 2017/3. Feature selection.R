# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Feature selection. Select the most variable genes.


library(scater)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(knitr)
library(gridExtra)

#-------------------------------------Load data----------------------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/3.Feature_Selection")

sce = readRDS("results/sce_norm.rds")


#---------------------------------Variance of the log-counts-------------------------------------------------------

#Modeling Gene variation
dec.sce = modelGeneVar(sce)

# Visualizing the fit:
fit.sce = metadata(dec.sce)
jpeg(file="./figures/3.Feature_Selection/Mean_log_vs_variance.jpeg", width=6, height=4, units="in", res=300)
plot(fit.sce$mean, fit.sce$var,
     xlab="Mean of log-expression",
     ylab="Variance of log-expression",col="#008287")
curve(fit.sce$trend(x), col="#F7766F", add=TRUE,
      lwd=2)
dev.off()

# Ordering by most interesting genes for inspection.

dec.sce.a = as_tibble(dec.sce[order(dec.sce$bio, decreasing=TRUE),][1:8,],rownames="rownames")
dec.sce.a

#------------------------------------Coefficient of variation----------------------------------------------------

dec.cv2 = modelGeneCV2(sce)
fit.cv2 = metadata(dec.cv2)

jpeg(file="./figures/3.Feature_Selection/Coeff_of_variation.jpeg", width=6, height=4, units="in", res=300)
plot(fit.cv2$mean, fit.cv2$cv2, log="xy",
     xlab="Mean of log-expression",
     ylab="Squared coefficient of variation",col="#008287")
curve(fit.cv2$trend(x), col="#F7766F", add=TRUE, lwd=2)
dev.off()

dec.cv2[order(dec.cv2$ratio, decreasing=TRUE),]

#----------------------------------------HVGs-----------------------------------------------------------------

# Selecting HVGs on the largest metrics
hvg.sce.var = getTopHVGs(dec.sce, n=1000)
str(hvg.sce.var)

# Select 10% of most variable genes for downstream analysis.
chosen = getTopHVGs(dec.sce,prop=0.1)
str(chosen)

lapply(chosen, write, "HVGs.txt", append=TRUE, ncolumns=1000)

# Subsetting to just the HVGs --------------------------------------------------

sce.hvg = sce[chosen,]
dim(sce.hvg)

#-----------------------------------Save---------------------------------------

saveRDS(sce.hvg,"results/sce_hvg.rds")
saveRDS(sce,"results/sce_feature_selection.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
