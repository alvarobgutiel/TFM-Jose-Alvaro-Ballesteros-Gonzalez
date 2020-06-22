# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# Mathys et al., 2019
# Feature selection. Select the most variable genes.


library(scater)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(knitr)

#-------------------------------------Load data----------------------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/3.Feature_Selection")

sce = readRDS("results/sce_norm.rds")


#---------------------------------Variance of the log-counts-------------------------------------------------------

#Modeling Gene variation
dec.sce = modelGeneVar(sce)

# Visualizing the fit:
fit.sce = metadata(dec.sce)
tiff(file="./figures/3.Feature_Selection/Mean_log_vs_variance.tiff", width=6, height=4, units="in", res=300)
plot(fit.sce$mean, fit.sce$var,
     xlab="Mean of log-expression",
     ylab="Variance of log-expression",col="#440154")
curve(fit.sce$trend(x), col="#B8DE29", add=TRUE,
      lwd=2)
dev.off()

# Ordering by most interesting genes for inspection.
dec.sce[order(dec.sce$bio, decreasing=TRUE),]

#------------------------------------Coefficient of variation----------------------------------------------------

dec.cv2 = modelGeneCV2(sce)
fit.cv2 = metadata(dec.cv2)

tiff(file="./figures/3.Feature_Selection/Coeff_of_variation.tiff", width=6, height=4, units="in", res=300)
plot(fit.cv2$mean, fit.cv2$cv2, log="xy",col="#440154")
curve(fit.cv2$trend(x), col="#B8DE29", add=TRUE, lwd=2)
dev.off()

dec.cv2[order(dec.cv2$ratio, decreasing=TRUE),]

#------------------------------------Quantifying technical noise---------------------------------------------

dec.pois = modelGeneVarByPoisson(sce)
dec.pois = dec.pois[order(dec.pois$bio, decreasing=TRUE),]
head(dec.pois)

tiff(file="./figures/3.Feature_Selection/Mean_of_log_express_vs_Variance.tiff", width=6, height=4, units="in", res=300)
plot(dec.pois$mean, dec.pois$total, pch=16, xlab="Mean of log-expression",
     ylab="Variance of log-expression",col="#440154")
curve(metadata(dec.pois)$trend(x), col="#B8DE29", add=TRUE)
dev.off()

#----------------------------------Accounting for blocking factors-------------------------------------------

dec.block = modelGeneVar(sce, block=sce$batch)
head(dec.block[order(dec.block$bio, decreasing=TRUE),1:6])

par(mfrow=c(1,2))
blocked.stats= dec.block$per.block
for (i in colnames(blocked.stats)) {
  current = blocked.stats[[i]]
  tiff(file=paste0("./figures/3.Feature_Selection/Mean-log-expr_vs_Vari_log_exprss_",i,".tiff"), width=6, height=4, units="in", res=300)
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit = metadata(current)
  points(curfit$mean, curfit$var, col="red", pch=16)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
  dev.off()
}


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
