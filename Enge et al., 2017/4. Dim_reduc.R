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
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/4.Dimentional_reduction")

sce = readRDS("results/sce_feature_selection.rds")
sce.hvg = readRDS("results/sce_hvg.rds")

names(colData(sce))[names(colData(sce))=="donor_organism.development_stage.ontology_label"] = "Patients"
names(colData(sce))[names(colData(sce))=="donor_organism.sex"] = "Sex"

names(colData(sce.hvg))[names(colData(sce.hvg))=="donor_organism.development_stage.ontology_label"] = "Patients"
names(colData(sce.hvg))[names(colData(sce.hvg))=="donor_organism.sex"] = "Sex"

colfunc = colorRampPalette(c("#00383A", "#99E3E6"))
# ----------------------------Principal component analysis-------------------------------------------------

#Run PCA
sce.hvg = runPCA(sce.hvg,exprs_values = "logcounts")


# Using the elbow point

percent.var = attr(reducedDim(sce.hvg), "percentVar")
chosen.elbow = PCAtools::findElbowPoint(percent.var)
chosen.elbow

#Plotting PCA
reducedDim(sce.hvg, "PCA") = reducedDim(sce.hvg, "PCA")[,1:chosen.elbow]

jpeg(file="./figures/4.Dimentional_reduction/PCA_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by= "Sex")+ggtitle("PCA plot by sex")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_manual(values = colfunc(2))
dev.off()

jpeg(file="./figures/4.Dimentional_reduction/PCA_by_age.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="Patients")+ggtitle("PCA plot by age")+xlab("PCA1")+ylab("PCA2")+
  theme(legend.text = element_text(size=7))
dev.off()

jpeg(file="./figures/4.Dimentional_reduction/PCA_by_sex_4_comp.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", ncomponents=4,colour_by="donor_organism.sex")+ggtitle("PCA plot by sex with 4 components")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction/PCA_by_condition_4_age.jpeg", width=10, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", ncomponents=4,colour_by="donor_organism.development_stage.ontology_label")+ggtitle("PCA plot by age with 4 components")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction/PCA_by_sex_and_age.jpeg", width=9, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="donor_organism.development_stage.ontology_label",shape_by="donor_organism.sex")+ggtitle("PCA plot by age and sex")
dev.off()


#---------------------------------------------tSNE----------------------------------------------------------

sce.hvg = runTSNE(sce.hvg, dimred="PCA")

jpeg(file="./figures/4.Dimentional_reduction/tSNE_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="Sex")
dev.off()

jpeg(file="./figures/4.Dimentional_reduction/tSNE_by_age.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="Patients")+
  theme(legend.text = element_text(size=7))
dev.off()

#Perplexity (5,20,80)
sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=5)
out5 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="donor_organism.sex") + ggtitle("perplexity = 5")

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=20)
out20 = plotReducedDim(sce.hvg, dimred="TSNE",
                        colour_by="donor_organism.sex") + ggtitle("perplexity = 20")

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=80)
out80 = plotReducedDim(sce.hvg, dimred="TSNE",
                        colour_by="donor_organism.sex") + ggtitle("perplexity = 80")

jpeg(file="./figures/4.Dimentional_reduction/tSNE_multiplot_out5_20_80_sex.jpeg", width=10, height=4, units="in", res=300)
multiplot(out5, out20, out80, cols=3)
dev.off()


#Perplexity (5,20,80)
sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=5)
out5 = plotReducedDim(sce.hvg, dimred="TSNE",
                      colour_by="donor_organism.development_stage.ontology_label") + ggtitle("perplexity = 5")

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=20)
out20 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="donor_organism.development_stage.ontology_label") + ggtitle("perplexity = 20")

sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity=80)
out80 = plotReducedDim(sce.hvg, dimred="TSNE",
                       colour_by="donor_organism.development_stage.ontology_label") + ggtitle("perplexity = 80")

jpeg(file="./figures/4.Dimentional_reduction/tSNE_multiplot_out5_20_80_age.jpeg", width=18, height=4, units="in", res=300)
multiplot(out5, out20, out80, cols=3)
dev.off()

#----------------------------------------------UMAP-----------------------------------------------------------

sce.hvg = runUMAP(sce.hvg, dimred="PCA")

jpeg(file="./figures/4.Dimentional_reduction/UMAP_sex.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="Sex")+
  theme(legend.text = element_text(size=7))
dev.off()

jpeg(file="./figures/4.Dimentional_reduction/UMAP_age.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="Patients")+
  theme(legend.text = element_text(size=7))
dev.off()


#-----------------------------------------------Save----------------------------------------------------------

saveRDS(sce.hvg,"results/dim_reduc.rds")

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
