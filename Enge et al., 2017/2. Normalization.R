# March 2010

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Normalization

library(SingleCellExperiment)
library(scater)
library(scran)

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/2.Normalization")

#-------------------Load data---------------------------------------------------

sce = readRDS("results/sce_qc.rds")

names(colData(sce))[names(colData(sce))=="donor_organism.development_stage.ontology_label"] = "Patients"
names(colData(sce))[names(colData(sce))=="donor_organism.sex"] = "Sex"
#------------------------------Normalization-------------------------------------

quick_clusters = quickCluster(sce, use.ranks = FALSE)

sce = computeSumFactors(sce,clusters = quick_clusters)
sce = logNormCounts(sce)

sce@colData$nCountnomr = Matrix::colSums(logcounts(sce))
metadata_4 = as.data.frame(sce@colData)

summary(sizeFactors(sce))

jpeg(file="./figures/2.Normalization/Size_factors_effect.jpeg", width=6, height=4, units="in", res=300)
hist(sizeFactors(sce))
dev.off()

colfunc = colorRampPalette(c("#F7766F", "#00BBC2"))

jpeg(file="./figures/2.Normalization/NCells_vs_NGenes_normalized_1.jpeg", width=8, height=6, units="in", res=300)
metadata_4 %>%
  ggplot(aes(x=Patients, y=nCountnomr, fill=Patients)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Counts per patient for normalized data")+
  scale_fill_manual(values = colfunc(8))
dev.off()

jpeg(file="./figures/2.Normalization/NCells_vs_NGenes_raw_1.jpeg", width=8, height=6, units="in", res=300)
metadata_4 %>%
  ggplot(aes(x=Patients, y=nCount, fill=Patients)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Counts per patient for raw data")+
  scale_fill_manual(values = colfunc(8))
dev.off()

#-----------------------------------------------------------------------------

saveRDS(sce,"results/sce_norm.rds")

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
