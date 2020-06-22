# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Dimentional reduction


library(scater)
library(scran)
library(SingleCellExperiment)
library(cowplot)
library(ggplot2)
library(PCAtools)


#-------------------------------------Load data----------------------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/4.Dimentional_reduction")
dir.create("./figures/4.Dimentional_reduction_20")

sce = readRDS("results/sce_feature_selection.rds")
sce.hvg = readRDS("results/sce_hvg.rds")

# ----------------------------Principal component analysis-------------------------------------------------

#Run PCA
sce.hvg = runPCA(sce.hvg,exprs_values = "logcounts")

# Using the elbow point

percent.var = attr(reducedDim(sce.hvg), "percentVar")
chosen.elbow = PCAtools::findElbowPoint(percent.var)
chosen.elbow

#Plotting PCA
reducedDim(sce.hvg, "PCA") = reducedDim(sce.hvg, "PCA")[,1:chosen.elbow]

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_patient_20.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="patient")+
  ggtitle("PCA plot by patients")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="condition")+
  ggtitle("PCA plot by condition")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_patient_4_comp.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", ncomponents=4,colour_by="patient")+
  ggtitle("PCA plot by patients with 4 components")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_condition_4_comp.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", ncomponents=4,colour_by="condition")+
  ggtitle("PCA plot by condition with 4 components")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_sex.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="sex")+
  ggtitle("PCA plot by sex")+
  scale_fill_manual(values= c("#E15635","#2D0A5A"))
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_patient_and_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="patient",shape_by="condition")+
  ggtitle("PCA plot by patient and condition")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_sex_and_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="sex",shape_by="condition")+
  ggtitle("PCA plot by sex and condition")+
  scale_fill_viridis_d(option="B")
dev.off()


#---------------------------------------------tSNE----------------------------------------------------------

sce.hvg = runTSNE(sce.hvg, dimred="PCA")

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_by_patient.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_by_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="condition")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_by_sex1.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="sex")+
  scale_fill_manual(values= c("#2D0A5A","#FBBC21"))
dev.off()

#Perplexity (5,20,80)
sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=5)
out5 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="patient") + ggtitle("perplexity = 5")+
  scale_fill_viridis_d(option="B")

out5_cond = plotReducedDim(sce.hvg, dimred="TSNE",
                           colour_by="condition") + ggtitle("perplexity = 5")+
  scale_fill_viridis_d(option="B")

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=20)
out20 = plotReducedDim(sce.hvg, dimred="TSNE",
                        colour_by="patient") + ggtitle("perplexity = 20")+
  scale_fill_viridis_d(option="B")

out20_cond = plotReducedDim(sce.hvg, dimred="TSNE",
                            colour_by="condition") + ggtitle("perplexity = 20")+
  scale_fill_viridis_d(option="B")

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=80)
out80 = plotReducedDim(sce.hvg, dimred="TSNE",
                        colour_by="patient") + ggtitle("perplexity = 80")+
  scale_fill_viridis_d(option="B")
out80_cond = plotReducedDim(sce.hvg, dimred="TSNE",
                            colour_by="condition") + ggtitle("perplexity = 80")+
  scale_fill_viridis_d(option="B")

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_multiplot_out5_20_80.jpeg", width=15, height=4, units="in", res=300)
multiplot(out5, out20, out80, cols=3)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_multiplot_out5_20_80_condition.jpeg", width=15, height=4, units="in", res=300)
multiplot(out5_cond, out20_cond, out80_cond, cols=3)
dev.off()

#----------------------------------------------UMAP-----------------------------------------------------------

sce.hvg = runUMAP(sce.hvg, dimred="PCA")

jpeg(file="./figures/4.Dimentional_reduction_20/UMAP_patient.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/UMAP_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="condition")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/UMAP_sex.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="sex")+
  scale_fill_manual(values= c("#2D0A5A","#FBBC21"))
dev.off()


#-----------------------------------------------Save----------------------------------------------------------

saveRDS(sce.hvg,"results/dim_reduc_20.rds")

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
