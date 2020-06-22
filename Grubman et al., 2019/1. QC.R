
# March 2020

# Single cell analysis - SingleCellExperiment
# GSE138852
#QC analysis

#------------------------------Library---------------------------------------------------

library(SingleCellExperiment)
library(scater)
library(scran)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cowplot)


#--------------------------Make directories---------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)

#dir.create("data")
#dir.create("results")
dir.create("figures")
#dir.create("scripts")

#----------------------------Loading data--------------------------------------------------

#Read counts matrix and convert it to a matrix object. Change rownames to correct names
counts_matrix = read.csv("./data/GSE138852_counts.csv")
rownames(counts_matrix) = counts_matrix[,1] #Asignar los nombres de la primera columna a las celdas
counts_matrix = counts_matrix[,-1] #Eliminamos la primera columna si contiene el nombre de las celdas
counts_matrix = as.matrix(counts_matrix) #Convertimos a matriz

#Set SingleCellExperiment object with matrix
sce = SingleCellExperiment(assays = list(counts = counts_matrix)) #SingleCellExperiment object with counts matrix as assay


#---------------------------Preparing data----------------------------------------------

#Metadata
cell_metadata = read.csv("./data/scRNA_metadata.tsv",sep="\t") #Import metadata if available
rownames(cell_metadata) = cell_metadata[,1] #Assign cell names to row names
cell_metadata = cell_metadata[,-1] #Remove first column if cell names is the first column

colData(sce) = DataFrame(cell_metadata) #Assign metadata to colData of sce object

#-------------------------------Condition--------------------------------------------------

num.cells = length(colnames(sce)) #Number of cells in sce

#Assign conditions to each cell (in this case, AD vs Control)
for (i in 1:num.cells){
  if (grepl("AD",sce$patient[i]) == TRUE){
    sce$condition[i] <- "AD"
  }else{
    sce$condition[i] <- "Control"
}}

#------------------------------------QC---------------------------------------------------

sce@colData$nCount = Matrix::colSums(counts(sce)) #Library size. Number of UMI counts taken across all cells
sce@colData$nFeatures = Matrix::colSums(counts(sce)>0) #Cells with any counts for each gene.
sce@colData$size_factors = scater::librarySizeFactors(sce)
sce = addPerCellQC (sce)
sce = addPerFeatureQC (sce)

#Mitochondrial genes
mito_genes = rownames(sce)[grep("^MT-",rownames(sce))] #Detect genes related with mitochondrial functions
sce@colData$percent_mito = Matrix::colSums(counts(sce)[mito_genes, ]) / sce@colData$nCount #Percentage of mitochondrial genes per cell

#Ribosome genes
ribo_genes = rownames(sce)[grep("^RP[SL]",rownames(sce))] #Detect genes related with ribosome functions
sce@colData$percent_ribo = Matrix::colSums(counts(sce)[ribo_genes, ]) / sce@colData$nCount #Percentage of ribosome genes

#Create metadata dataframe
metadata_2 = as.data.frame(sce@colData)

#Plots of QC before filtering

dir.create("./figures/1.QC")
dir.create("./figures/1.QC/Before_filtering")

jpeg(file="./figures/1.QC/Before_filtering/density_plot_nFeatures.jpeg", width=6, height=4, units="in", res=300)
plot(density(sce$nFeatures),main="Density - total counts")
abline(v=200)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nFeatures_by_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "nFeatures",x = "patient",colour_by = "patient")+
  ggtitle("nFeatures per patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_by_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "nCount",x = "patient",colour_by = "patient")+
  ggtitle("nCount per patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_by_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "percent_mito",x = "patient",colour_by = "patient")+
  ggtitle("Percentage of mitochondrial genes per patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_vs._nFeatures.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "patient")+
  ggtitle("nCount vs. nFeatures")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_vs._nFeatures.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "patient")+
  ggtitle("Percentage of mitochondrial genes vs. nFeatures")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/NCells_vs_NGenes_filtered.jpeg", width=8, height=6, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=patient, y=nCount, fill=patient)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by patient")+
  scale_fill_viridis_d(option="B")
dev.off()

#Plot number of genes by patient
jpeg(file="./figures/1.QC/Before_filtering/NCells_filtered.jpeg", width=8, height=6, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=patient, fill=patient)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells by patient")+
  scale_fill_viridis_d(option="B")
dev.off()

# Visualize the number UMIs/transcripts per cell
jpeg(file="./figures/1.QC/Before_filtering/NUMIs_filtered.jpeg", width=8, height=6, units="in", res=300)
metadata_2 %>%
  ggplot(aes(color=patient, x=nCount, fill= patient))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by patient")+
  scale_fill_viridis_d(option="B")
dev.off()

#-------------------------------------Filtering--------------------------------------------

#Previous number of genes/cells. Data was pre-filtered.
dim(sce)

threshold = 200 #Threshold for nFeatures (number of genes detected in each cell)

selected_c =  colnames(sce)[sce$nFeatures > threshold]
selected_f = rownames(sce)[ Matrix::rowSums(counts(sce)) > 3]
sce.filt = sce[selected_f , selected_c]

#Check filtered genes/cells
dim(sce.filt)

high.det.v3 = sce.filt$nFeatures > 4100
# remove these cells
sce.filt = sce.filt[ , (!high.det.v3)]
# check number of cells
ncol(sce.filt)

#Compute the relative expression of each gene per cell
rel_expression = t(t(counts(sce.filt)) / Matrix::colSums(counts(sce.filt))) * 100
most_expressed = sort(Matrix::rowSums(rel_expression),T)[20:1] / ncol(sce.filt)

jpeg(file="./figures/1.QC/most_expressed.jpeg", width=6, height=4, units="in", res=300)
boxplot(as.matrix(t(rel_expression[names(most_expressed),])),cex=.1, las=1, xlab="% total count per cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()

#Select Mitoratio under 0.20%
selected_mito = sce.filt$percent_mito < 0.20
# and subset the object to only keep those cells
sce.filt = sce.filt[, selected_mito]
dim(sce.filt)

#--------------------------Plot after filtering------------------------------------------------------------

dir.create("./figures/1.QC/After_filtering")

#Create metadata dataframe
metadata_3 = as.data.frame(sce.filt@colData)

jpeg(file="./figures/1.QC/After_filtering/density_plot_nFeatures_AF.jpeg", width=6, height=4, units="in", res=300)
plot(density(sce.filt$nFeatures),main="Density - total counts")
abline(v=200)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nFeatures_AF.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,y = "nFeatures",x = "patient",colour_by = "patient")+
  ggtitle("nFeatures after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCounts_AF.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,y = "nCount",x = "patient",colour_by = "patient")+
  ggtitle("nCounts after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_AF.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,y = "percent_mito",x = "patient",colour_by = "patient")+
  ggtitle("Percentage of mitochondrial genes after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCounts_vs._nFeatures_AF.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,x = "nCount",y = "nFeatures",colour_by = "patient")+
  ggtitle("nCounts vs. nFeatures after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_vs_nFeatures_AF.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,x = "percent_mito",y = "nFeatures",colour_by = "patient")+
  ggtitle("Percentage of mitochondrial genes vs. nFeatures after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/NCells_vs_NGenes_filtered_AF.jpeg", width=8, height=6, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=patient, y=log10(nCount), fill=patient)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by patient after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

#Plot number of genes by patient
jpeg(file="./figures/1.QC/After_filtering/NCells_filtered_AF.jpeg", width=8, height=6, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=patient, fill=patient)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells by patient after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

# Visualize the number UMIs/transcripts per cell
jpeg(file="./figures/1.QC/After_filtering/NUMIs_filtered_AF.jpeg", width=8, height=6, units="in", res=300)
metadata_3 %>%
  ggplot(aes(color=patient, x=nCount, fill= patient))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by patient after filtering")+
  scale_fill_viridis_d(option="B")
dev.off()

#-------------------------------Cell cycle--------------------------------------------------

hm.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
symb = mapIds(org.Hs.eg.db, keys=rownames(sce.filt), keytype="SYMBOL", column="ENSEMBL")
rowData(sce.filt)$SYMBOL = rownames(sce.filt)
rowData(sce.filt)$ENSEMBL = symb
assignments = cyclone(sce.filt, hm.pairs, gene.names=rowData(sce.filt)$ENSEMBL)
sce.filt$CellCycle = assignments$phases

#Scores
sce.filt$G1.score = assignments$scores$G1
sce.filt$G2M.score = assignments$scores$G2M
sce.filt$S.score = assignments$scores$S

#Plot scores
dir.create("./figures/1.QC/Cell_cycle")

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score_vs_G1_score.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "G1.score",colour_by = "patient")+
  ggtitle("G2M score vs. G1 score")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "patient",colour_by = "patient")+
  ggtitle("G2M score by patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G1_score.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G1.score",x = "patient",colour_by = "patient")+
  ggtitle("G1 score by patient")+
  scale_fill_viridis_d(option="B")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/S_score.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "S.score",x = "patient",colour_by = "patient")+
  ggtitle("S score by patient")+
  scale_fill_viridis_d(option="B")
dev.off()


#-----------------------------Save results----------------------------------------------

saveRDS(sce.filt,"results/sce_qc.rds")

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

