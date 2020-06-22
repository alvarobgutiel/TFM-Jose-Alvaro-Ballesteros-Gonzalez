# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# Mathys et al., 2019
# Clustering

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(umap)
library(pheatmap)
library(igraph)
library(cowplot)
library(SingleR)
library(BRETIGEA)
library(reshape2)
library(ggpubr)
library(viridis)
library(Matrix)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/5.Clustering_20")

sce.hvg = readRDS("results/dim_reduc_20.rds")


#-------------------------------------------SNN clustering----------------------------------------------------

g = buildSNNGraph(sce.hvg,k=30,use.dimred="PCA")

cluster = igraph::cluster_walktrap(g)$membership
sce.hvg$cluster = factor(cluster)
table(sce.hvg$cluster)

#tSNE by clusters
jpeg(file="./figures/5.Clustering_20/tSNE_by_clusters.jpeg", width=10, height=6, units="in", res=300)
plotTSNE(sce.hvg,colour_by="cluster", text_by="cluster")+
  scale_fill_viridis_d(direction = -1)
dev.off()

#PCA by clusters
jpeg(file="./figures/5.Clustering_20/PCA_by_clusters.jpeg", width=10, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="cluster", text_by="cluster")+
  scale_fill_viridis_d(direction = -1)
dev.off()

#UMAP by clusters
jpeg(file="./figures/5.Clustering_20/UMAP_by_clusters.jpeg", width=10, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="cluster", text_by="cluster")+
  scale_fill_viridis_d(direction = -1)
dev.off()

#Ratio of observed to expected edge
ratio = clusterModularity(g, cluster, as.ratio=TRUE)

jpeg(file="./figures/5.Clustering_20/Pheatmap1.jpeg", width=10, height=6, units="in", res=300)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=viridis(100))
dev.off()

#Bootstrap to compute the co-assignment probability
ass.prob = bootstrapCluster(sce.hvg, FUN=function(x) {
  g = buildSNNGraph(x, use.dimred="PCA")
  igraph::cluster_walktrap(g)$membership
}, clusters=sce.hvg$cluster)

jpeg(file="./figures/5.Clustering_20/Pheatmap2.jpeg", width=10, height=6, units="in", res=300)
pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=viridis(100))
dev.off()

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
 
 jpeg(file="./figures/5.Clustering_20/PCA_HC_euclidean_correlation.jpeg", width=12, height=6, units="in", res=300)
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
#-------------------------------------------BRETIGEA--------------------------------------------------------------------

sce.matrix = as.matrix(logcounts(sce.hvg))

ct_res = brainCells(sce.matrix, nMarker = 50)
type = colnames(ct_res)[apply(ct_res,1,which.max)]
sce.hvg$cell.type = type

sce.hvg$cell.type = factor(sce.hvg$cell.type, levels = c("ast","end","mic","neu","oli","opc"),labels = c("Astrocytes", "Endothelial Cells", "Microglia",
                               "Neurons", "Oligodendrocytes", "Oligodendrocyte PC"))

ct_res_df = as.data.frame(ct_res)
ct_res_df$type = sce.hvg$pathology.group[rownames(colData(sce.hvg)) == rownames(ct_res_df)]
ct_res_df = melt(ct_res_df)

jpeg(file="./figures/5.Clustering_20/BRETIGEA_score.jpeg", width=6, height=4, units="in", res=300)
ggplot(ct_res_df, aes(x = type, y = value, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    facet_wrap(~variable) +
    xlab("") + ylab("BRETIGEA score") +
    guides(fill = FALSE)+
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_viridis_d()
dev.off()

#---------------------------------------Plots with classes------------------------------------------------------------

jpeg(file="./figures/5.Clustering_20/PCA_by_cell_type.jpeg", width=12, height=8, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="cell.type")+
  ggtitle("PCA plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/5.Clustering_20/tSNE_by_cell_type.jpeg", width=12, height=8, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="cell.type")+
  ggtitle("tSNE plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/5.Clustering_20/UMAP_by_cell_type.jpeg", width=12, height=8, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP", colour_by="cell.type")+
  ggtitle("UMAP plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(direction = -1)
dev.off()

#--------------------------------------------------------------------------------------------------------

jpeg(file="./figures/5.Clustering_20/Expression_Genes.jpeg", width=8, height=10, units="in", res=300)
plotExpression(sce.hvg, features=c("APOE","GRID2","NRXN1","RBFOX1","PLP1","NRGN","GAD1","ERBB4"),
               x=I(factor(sce.hvg$cell.type)), colour_by=I(factor(sce.hvg$cell.type)))+
  scale_fill_viridis_d(direction=-1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=6))
dev.off()


#------------------------------------------------------Save------------------------------------------------------------

saveRDS(sce.hvg,"results/clustering_20.rds")

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
