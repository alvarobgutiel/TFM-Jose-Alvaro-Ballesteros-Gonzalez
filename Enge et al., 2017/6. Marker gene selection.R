# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Marker gene detection

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(rafalib)
library(MAST)
library(tibble)
library(dplyr)

#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/6.Marker Gene Selection")

sce.hvg = readRDS("results/clustering.rds")

#-------------------------------------------------------------------------------------------------------------

#Compute differential expression
markers_genes = findMarkers( x = sce.hvg,
                              groups = sce.hvg$cluster,
                              lfc=.5,
                              pval.type = "all",
                              direction = "up")

#Colect the top 25 genes for each cluster and put them into a single table
top25 = lapply( names(markers_genes), function(x) { temp = markers_genes[[x]][1:25, 1:2] ; temp$gene = rownames(markers_genes[[x]])[1:25] ; temp$cluster = x ; return(temp) } )
top25 = as_tibble(do.call(rbind, top25))
write.csv(top25, file = "top25.csv")


#25 gene expression of each cluster

for(i in unique(top25$cluster)){
  tiff(file=paste0("./figures/6.Marker Gene Selection/top25_cluster",i,".tiff"), width=6, height=4, units="in", res=300)
  barplot( sort( setNames(-log10(top25$p.value), top25$gene) [top25$cluster == i], F),
           horiz = T,las=1 ,main=paste0(i," vs. rest"),border = "white", yaxs="i",xlab="-log10FC" )
  abline(v=c(0,-log10(0.05)),lty=c(1,2))
  dev.off()
}

colfunc = colorRampPalette(c("#F7766F", "#005c60"))

#HeatMap
as_tibble(top25) %>% group_by(cluster)  %>% top_n(-5, p.value) -> top5
tiff(file="./figures/6.Marker Gene Selection/Heatmap.tiff", width=6, height=6, units="in", res=300)
scater::plotHeatmap(sce.hvg[,order(sce.hvg$cluster)], features = unique(top5$gene) ,
                    center=T , zlim = c(-3,3) ,
                    colour_columns_by= "cluster",
                    show_colnames=F , cluster_cols=F,
                    fontsize_row=6,
                    color=colfunc(90)
)
dev.off()

#----------------------------------------------------------------------------------------------------------------------

#Compute differential expression
markers_genes_cell_type = findMarkers( x = sce.hvg,
                             groups = sce.hvg$cell.type,
                             lfc=.5,
                             pval.type = "all",
                             direction = "up")

#Colect the top 25 genes for each cluster and put them into a single table
top25_cell_type = lapply( names(markers_genes_cell_type), function(x) { temp = markers_genes_cell_type[[x]][1:25, 1:2] ; temp$gene = rownames(markers_genes_cell_type[[x]])[1:25] ; temp$cell.type = x ; return(temp) } )
top25_cell_type = as_tibble(do.call(rbind, top25_cell_type))
write.csv(top25_cell_type, file = "top25_celltype.csv")


#25 gene expression of each cluster
colfunc = colorRampPalette(c("#F7766F", "#005c60"))

#HeatMap
as_tibble(top25_cell_type) %>% group_by(cell.type)  %>% top_n(-5, p.value) -> top5_cell_type
jpeg(file="./figures/6.Marker Gene Selection/Heatmap_cell_type_1.jpeg", width=6, height=6, units="in", res=300)
scater::plotHeatmap(sce.hvg[,order(sce.hvg$cell.type)], features = unique(top5_cell_type$gene) ,
                    center=T , zlim = c(-3,3) ,
                    colour_columns_by= "cell.type",
                    show_colnames=F , cluster_cols=F,
                    fontsize_row=6,
                    color=colfunc(90),
                    main="Heatmap per cell type"
)
dev.off()
#----------------------------Save-------------------------------------------------------------------------

saveRDS(sce.hvg,"results/mark_gene_select.rds")

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------



