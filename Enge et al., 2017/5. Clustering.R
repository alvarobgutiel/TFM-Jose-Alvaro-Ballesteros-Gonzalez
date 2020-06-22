# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Clustering

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(umap)
library(pheatmap)
library(igraph)
library(cowplot)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/5.Clustering")

sce.hvg = readRDS("results/dim_reduc.rds")


#-------------------------------------------SNN clustering----------------------------------------------------

g = buildSNNGraph(sce.hvg,k=20,use.dimred="PCA",assay.type="RNA")

cluster = igraph::cluster_walktrap(g)$membership
sce.hvg$cluster = factor(cluster)
table(sce.hvg$cluster)

#tSNE by clusters
jpeg(file="./figures/5.Clustering/tSNE_by_clusters.jpeg", width=6, height=4, units="in", res=300)
plotTSNE(sce.hvg,colour_by="cluster", text_by="cluster")
dev.off()

#PCA by clusters
jpeg(file="./figures/5.Clustering/PCA_by_clusters.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="cluster", text_by="cluster")
dev.off()

#UMAP by clusters
jpeg(file="./figures/5.Clustering/UMAP_by_clusters.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="cluster", text_by="cluster")
dev.off()

#Ratio of observed to expected edge
ratio = clusterModularity(g, cluster, as.ratio=TRUE)

jpeg(file="./figures/5.Clustering/Pheatmap1.jpeg", width=6, height=4, units="in", res=300)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=rev(heat.colors(100)))
dev.off()

#Bootstrap to compute the co-assignment probability
ass.prob = bootstrapCluster(sce.hvg, FUN=function(x) {
  g = buildSNNGraph(x, use.dimred="PCA")
  igraph::cluster_walktrap(g)$membership
}, clusters=sce.hvg$cluster)

jpeg(file="./figures/5.Clustering/Pheatmap2.jpeg", width=6, height=4, units="in", res=300)
pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=colorRampPalette(c("white", "blue"))(100))
dev.off()

#Subclusters
subout = quickSubCluster(sce.hvg, sce.hvg$cluster)

#----------------------------------------Correlations-------------------------------------------------------

d = dist(reducedDims(sce.hvg)$PCA,
           method="euclidean")
#Compute sample correlations
sample_cor = cor( Matrix::t(reducedDims(sce.hvg)$PCA) )
#Transform the scale from correlations
sample_cor = (1 - sample_cor) / 2
#Convert it to a distance object
d2 = as.dist(sample_cor)

#euclidean
h_euclidean = hclust(d, method="ward.D2")
#correlation
h_correlation = hclust(d2, method="ward.D2")

#euclidean distance
sce.hvg$hc_euclidean_5 = factor( cutree(h_euclidean,k = 5) )
sce.hvg$hc_euclidean_10 = factor( cutree(h_euclidean,k = 10) )
sce.hvg$hc_euclidean_15 = factor( cutree(h_euclidean,k = 15) )
#correlation distance
sce.hvg$hc_corelation_5 = factor( cutree(h_correlation,k = 5) )
sce.hvg$hc_corelation_10 = factor( cutree(h_correlation,k = 10) )
sce.hvg$hc_corelation_15 = factor( cutree(h_correlation,k = 15) )

jpeg(file="./figures/5.Clustering/PCA_HC_euclidean_correlation.jpeg", width=12, height=6, units="in", res=300)
plot_grid(ncol = 3,
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_euclidean_5")+
            ggplot2::ggtitle(label ="HC_euclidean_5"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_euclidean_10")+
            ggplot2::ggtitle(label ="HC_euclidean_10"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_euclidean_15")+
            ggplot2::ggtitle(label ="HC_euclidean_15"),

          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_corelation_5")+
            ggplot2::ggtitle(label ="HC_correlation_5"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_corelation_10")+
            ggplot2::ggtitle(label ="HC_correlation_10"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_corelation_15")+
            ggplot2::ggtitle(label ="HC_correlation_15")
)
dev.off()

#-----------------------------------------SingleR---------------------------------------------------------

library(SingleR)
ref = BlueprintEncodeData()
pred = SingleR(test=sce.hvg, ref=ref, labels=ref$label.main)
table(pred$labels)
sce.hvg$cell.type = pred$labels

subset_celltypes1 = sce.hvg[,sce.hvg$cell.type == "B-cells"]
subset_celltypes2 = sce.hvg[,sce.hvg$cell.type == "Epithelial cells"]
subset_celltypes3 = sce.hvg[,sce.hvg$cell.type == "Fibroblast"]
subset_celltypes4 = sce.hvg[,sce.hvg$cell.type == "HSC"]
subset_celltypes5 = sce.hvg[,sce.hvg$cell.type == "Monocytes"]
subset_celltypes6 = sce.hvg[,sce.hvg$cell.type == "Neurons"]

subset_celltypes = cbind(subset_celltypes1,subset_celltypes2,subset_celltypes3,subset_celltypes4,subset_celltypes5,subset_celltypes6)
names(colData(subset_celltypes))[names(colData(subset_celltypes))=="cell.type"] = "Cell type"

#tSNE by clusters
jpeg(file="./figures/5.Clustering/tSNE_by_cell_type.jpeg", width=6, height=4, units="in", res=300)
plotTSNE(subset_celltypes,colour_by="Cell type",text_by="Cell type",text_size = 2.5)+
  theme(legend.text = element_text(size=10))
dev.off()

#PCA by clusters
jpeg(file="./figures/5.Clustering/PCA_by_cell_type.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(subset_celltypes, dimred="PCA",colour_by="Cell type",text_by="Cell type",text_size = 2.5)+
  theme(legend.text = element_text(size=10),)
dev.off()

#UMAP by clusters
jpeg(file="./figures/5.Clustering/UMAP_by_cell_type.jpeg", width=6, height=4, units="in", res=300)
plotReducedDim(subset_celltypes, dimred="UMAP",colour_by="Cell type",text_by="Cell type",text_size = 1.8)+
  theme(legend.text = element_text(size=10))
dev.off()



#--------------------------------------------Save------------------------------------------------------------

saveRDS(sce.hvg,"results/clustering.rds")

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
