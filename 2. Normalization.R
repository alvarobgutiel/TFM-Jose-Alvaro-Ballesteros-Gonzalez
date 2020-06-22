# March 2010

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Normalization

library(SingleCellExperiment)
library(scater)
library(scran)

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/2.Normalization")

#-------------------Load data---------------------------------------------------

sce = readRDS("results/sce_qc.rds")

# ------------------------------------------Library size factors ---------------------------------------------------------

lib.sce = librarySizeFactors(sce)

# Examine distribution of size factors
summary(lib.sce)

png(file="./figures/2.Normalization/Log10[Size_factor].png", width=6, height=4, units="in", res=300)
hist(log10(lib.sce), xlab="Log10[Size factor]", col="grey80")
dev.off()

ls.sce = colSums(counts(sce))

png(file="./figures/2.Normalization/Library_size_vs._Size_factor.png", width=6, height=4, units="in", res=300)
plot(ls.sce, lib.sce, log="xy", xlab="Library size", ylab="Size factor")
dev.off()

#------------------------------Normalization-------------------------------------
memory.limit(9999999999)
quick_clusters = quickCluster(sce, use.ranks = FALSE)

sce = computeSumFactors(sce,clusters = quick_clusters)
sce = logNormCounts(sce)

summary(sizeFactors(sce))

png(file="./figures/2.Normalization/Size_factors_effect.png", width=6, height=4, units="in", res=300)
hist(sizeFactors(sce))
dev.off()

----------------------------------plot norm data--------------------------------------------

sce@colData$nCountnomr = Matrix::colSums(logcounts(sce))
metadata_4 = as.data.frame(sce@colData)

png(file="./figures/2.Normalization/NCells_vs_NGenes_normalized.png", width=8, height=6, units="in", res=300)
metadata_4 %>%
  ggplot(aes(x=patient, y=nCountnomr, fill=patient)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by patient after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

#-----------------------------------------------------------------------------

saveRDS(sce,"results/sce_norm.rds")

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
