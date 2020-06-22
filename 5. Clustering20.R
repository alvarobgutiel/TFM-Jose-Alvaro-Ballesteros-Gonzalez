# March 2010

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
library(BRETIGEA)
library(reshape2)
library(ggpubr)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/5.Clustering_20")

sce.hvg = readRDS("results/dim_reduc_20.rds")


#-------------------------------------------SNN clustering----------------------------------------------------

g = buildSNNGraph(sce.hvg,k=20,use.dimred="PCA",assay.type="RNA")

cluster = igraph::cluster_walktrap(g)$membership
sce.hvg$cluster = factor(cluster)
table(sce.hvg$cluster)

#tSNE by clusters
jpeg(file="./figures/5.Clustering_20/tSNE_by_clusters.jpeg", width=8, height=6, units="in", res=300)
plotTSNE(sce.hvg,colour_by="cluster", text_by="cluster")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/tSNE_by_clusters_by_condition.jpeg", width=8, height=6, units="in", res=300)
plotTSNE(sce.hvg,colour_by="cluster", text_by="cluster",shape_by="condition")+
  scale_fill_viridis_d(option="B")
dev.off()

#PCA by clusters
jpeg(file="./figures/5.Clustering_20/PCA_by_clusters.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="cluster", text_by="cluster")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/PCA_by_clusters_by_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA",colour_by="cluster", text_by="cluster",shape_by="condition")+
  scale_fill_viridis_d(option="B")
dev.off()

#UMAP by clusters
jpeg(file="./figures/5.Clustering_20/UMAP_by_clusters.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="cluster", text_by="cluster")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/UMAP_by_clusters_by_conditon.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP",colour_by="cluster", text_by="cluster",shape_by="condition")+
  scale_fill_viridis_d(option="B")
dev.off()

#Ratio of observed to expected edge
ratio = clusterModularity(g, cluster, as.ratio=TRUE)

jpeg(file="./figures/5.Clustering_20/Pheatmap1.jpeg", width=6, height=4, units="in", res=300)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=viridis::inferno(100))
dev.off()

#Bootstrap to compute the co-assignment probability
ass.prob = bootstrapCluster(sce.hvg, FUN=function(x) {
  g = buildSNNGraph(x, use.dimred="PCA")
  igraph::cluster_walktrap(g)$membership
}, clusters=sce.hvg$cluster)

jpeg(file="./figures/5.Clustering_20/Pheatmap2.jpeg", width=6, height=4, units="in", res=300)
pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=viridis::inferno(100))
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

jpeg(file="./figures/5.Clustering_20/PCA_HC_euclidean_correlation.jpeg", width=12, height=6, units="in", res=300)
plot_grid(ncol = 3,
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_euclidean_5")+
            ggplot2::ggtitle(label ="HC_euclidean_5")+
            scale_fill_viridis_d(option="B"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_euclidean_10")+
            ggplot2::ggtitle(label ="HC_euclidean_10")+
            scale_fill_viridis_d(option="B"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_euclidean_15")+
            ggplot2::ggtitle(label ="HC_euclidean_15")+
            scale_fill_viridis_d(option="B"),

          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_corelation_5")+
            ggplot2::ggtitle(label ="HC_correlation_5")+
            scale_fill_viridis_d(option="B"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_corelation_10")+
            ggplot2::ggtitle(label ="HC_correlation_10")+
            scale_fill_viridis_d(option="B"),
          plotReducedDim(sce.hvg,dimred = "PCA",colour_by = "hc_corelation_15")+
            ggplot2::ggtitle(label ="HC_correlation_15")+
            scale_fill_viridis_d(option="B")
)
dev.off()

#-------------------------------------------BRETIGEA--------------------------------------------------------------------

sce.data.frame = as.data.frame(counts(sce.hvg))

ct_res = brainCells(sce.data.frame, nMarker = 50)
type = colnames(ct_res)[apply(ct_res,1,which.max)]
sce.hvg$cell.type = type

sce.hvg$cell.type = factor(sce.hvg$cell.type, levels = c("ast","end","mic","neu","oli","opc"),labels = c("Astrocytes", "Endothelial Cells", "Microglia",
                               "Neurons", "Oligodendrocytes", "Oligodendrocyte PC"))

ct_res_df = as.data.frame(ct_res)
ct_res_df = cbind(ct_res_df, colsplit(rownames(ct_res_df), pattern = "_", names = c("sample", "type")))
ct_res_df = melt(ct_res_df)

ct_res_df$variable = factor(ct_res_df$variable, levels = c("ast","end","mic","neu","oli","opc"),labels = c("Astrocytes", "Endothelial Cells", "Microglia",
                                                                                                           "Neurons", "Oligodendrocytes", "Oligodendrocyte PC"))

jpeg(file="./figures/5.Clustering_20/BRETIGEA_score.jpeg", width=8, height=6, units="in", res=300)
ggplot(ct_res_df, aes(x = type, y = value, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    facet_wrap(~variable) +
    xlab("") + ylab("BRETIGEA score") +
    guides(fill = FALSE)+
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_viridis_d(option="B")
dev.off()

#---------------------------------------Plots with classes------------------------------------------------------------

jpeg(file="./figures/5.Clustering_20/PCA_by_cell_type_1.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="cell.type")+
  ggtitle("PCA plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/tSNE_by_cell_type.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="cell.type")+
  ggtitle("tSNE plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/UMAP_by_cell_type.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP", colour_by="cell.type")+
  ggtitle("UMAP plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()


-------------------------------------------------------------------------------------------------------------------

jpeg(file="./figures/5.Clustering_20/PCA_by_cell_type_by_conditon.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="cell.type",shape_by="condition")+
  ggtitle("PCA plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/tSNE_by_cell_type_by_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="TSNE", colour_by="cell.type",shape_by="condition")+
  ggtitle("tSNE plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/5.Clustering_20/UMAP_by_cell_type_by_condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce.hvg, dimred="UMAP", colour_by="cell.type",shape_by="condition")+
  ggtitle("UMAP plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  scale_fill_viridis_d(option="B")
dev.off()

--------------------------------------------------------------------------------------------------------------------
subset_ad =  sce.hvg[,sce.hvg$condition == "AD"]
subset_cont =  sce.hvg[,sce.hvg$condition == "Control"]

PCA_ad = reducedDim(subset_ad)
tsne_ad = reducedDim(subset_ad,"TSNE")
umap_ad = reducedDim(subset_ad,"UMAP")
ad_types = subset_ad$cell.type

total_ad_pca = data.frame(PCA_ad,ad_types)
total_ad_tsne = data.frame(tsne_ad,ad_types)
total_ad_umap = data.frame(umap_ad,ad_types)

cont_PCA = data.frame(reducedDim(subset_cont))


jpeg(file="./figures/5.Clustering_20/PCA_by_cell_type_grey.jpeg", width=8, height=6, units="in", res=300)
ggplot(total_ad_pca,aes(x=PC1,y=PC2)) +
  geom_point(data=cont_PCA,aes(PC1,PC2),color="#DBDAD6")+
  geom_point(aes(colour = ad_types))+
  scale_colour_viridis_d()+
  ggtitle("PCA plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

jpeg(file="./figures/5.Clustering_20/tSNE_by_cell_type_grey.jpeg", width=8, height=6, units="in", res=300)
ggplot(total_ad_tsne,aes(x=X1,y=X2)) +
  geom_point(data=cont_PCA,aes(PC1,PC2),color="#DBDAD6",size=0.7)+
  geom_point(aes(colour = ad_types),size=0.7)+
  scale_colour_viridis_d()+
  ggtitle("tSNE plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

jpeg(file="./figures/5.Clustering_20/UMAP_by_cell_type_grey.jpeg", width=8, height=6, units="in", res=300)
ggplot(total_ad_umap,aes(x=X1,y=X2)) +
  geom_point(data=cont_PCA,aes(PC1,PC2),color="#DBDAD6",size=0.7)+
  geom_point(aes(colour = ad_types),size=0.7)+
  scale_colour_viridis_d()+
  ggtitle("UMAP plot by cell type")+xlab("PCA1")+ylab("PCA2")+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

#---------------------------------------------------------------------------------------------------------------------

jpeg(file="./figures/5.Clustering_20/Expression_Genes.jpeg", width=8, height=10, units="in", res=300)
plotExpression(sce.hvg, features=c("LINGO1","APOE","NEAT1","GRID2","GABRB1","NRXN1","RBFOX1","ERBB4"),
               x=I(factor(sce.hvg$cell.type)), colour_by=I(factor(sce.hvg$cell.type)))+
  scale_fill_viridis_d(option="B")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=6))
dev.off()

#------------------------------------------------------Save------------------------------------------------------------

saveRDS(sce.hvg,"results/clustering_20.rds")

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------