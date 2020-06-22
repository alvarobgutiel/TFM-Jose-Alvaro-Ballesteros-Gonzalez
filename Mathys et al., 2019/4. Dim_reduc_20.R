# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# Mathys et al., 2019
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
dir.create("./figures/4.Dimentional_reduction_20")

sce = readRDS("results/sce_feature_selection_20.rds")
sce.hvg = readRDS("results/sce_hvg_20.rds")

# ----------------------------Principal component analysis-------------------------------------------------

#Run PCA
sce.hvg = runPCA(sce.hvg,exprs_values = "logcounts")

# Using the elbow point

percent.var = attr(reducedDim(sce.hvg), "percentVar")
chosen.elbow = PCAtools::findElbowPoint(percent.var)
chosen.elbow

#Plotting PCA
reducedDim(sce.hvg, "PCA") = reducedDim(sce.hvg, "PCA")[,1:chosen.elbow]

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_pathology_group.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="pathology.group")+
  ggtitle("PCA plot by pathology group")+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_msex.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="msex")+ggtitle("PCA plot by sex")+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_path_group_4_comp.jpeg", width=10, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", ncomponents=4,colour_by="pathology.group")+
  ggtitle("PCA plot by pathology group with 4 components")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_msex_4_comp.jpeg", width=10, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", ncomponents=4,colour_by="msex")+
  ggtitle("PCA plot by sex with 4 components")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_unros.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="un_ros")+
  ggtitle("PCA plot by patient")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_pathology_group_and_sex.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="pathology.group",shape_by="msex")+
  ggtitle("PCA plot by pathology group and sex")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/PCA_by_patient_and_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="un_ros",shape_by="pathology.group")+
  ggtitle("PCA plot by patient and pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()


#---------------------------------------------tSNE----------------------------------------------------------

sce.hvg = runTSNE(sce.hvg, dimred="PCA")

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_by_pathology_group.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="pathology.group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_by_sex_2.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="msex")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_by_patient.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="un_ros")+
  scale_fill_viridis_d(direction = -1)
dev.off()

#Perplexity (5,20,80)
sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=5)
out5 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="pathology.group") + ggtitle("perplexity = 5")+
  scale_fill_viridis_d(direction = -1)

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=20)
out20 = plotReducedDim(sce.hvg, dimred="TSNE",
                        colour_by="pathology.group") + ggtitle("perplexity = 20")+
  scale_fill_viridis_d(direction = -1)

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=80)
out80 = plotReducedDim(sce.hvg, dimred="TSNE",
                        colour_by="pathology.group") + ggtitle("perplexity = 80")+
  scale_fill_viridis_d(direction = -1)

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_multiplot_out5_20_80.jpeg", width=14, height=4, units="in", res=300)
multiplot(out5, out20, out80, cols=3)
dev.off()


#Perplexity (5,20,80)
sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=5)
out5 = plotReducedDim(sce.hvg, dimred="TSNE",
                      colour_by="msex") + ggtitle("perplexity = 5")+
  scale_fill_viridis_d(direction = -1)

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=20)
out20 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="msex") + ggtitle("perplexity = 20")+
  scale_fill_viridis_d(direction = -1)

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=80)
out80 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="msex") + ggtitle("perplexity = 80")+
  scale_fill_viridis_d(direction = -1)

jpeg(file="./figures/4.Dimentional_reduction_20/tSNE_multiplot_out5_20_80_sex.jpeg", width=14, height=4, units="in", res=300)
multiplot(out5, out20, out80, cols=3)
dev.off()

#----------------------------------------------UMAP-----------------------------------------------------------

sce.hvg = runUMAP(sce.hvg, dimred="PCA")

jpeg(file="./figures/4.Dimentional_reduction_20/UMAP_pathology_group.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="pathology.group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/UMAP_un_ros.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="un_ros")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/4.Dimentional_reduction_20/UMAP_sex.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="msex")+
  scale_fill_viridis_d(direction = -1)
dev.off()


#-----------------------------------------------Save----------------------------------------------------------

saveRDS(sce.hvg,"results/dim_reduc_20.rds")

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
