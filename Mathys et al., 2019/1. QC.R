
# March 2010

# Single cell analysis - SingleCellExperiment
# GSE138852
#QC analysis

#------------------------------------------------------------------------------------------

#Metadata
#nCount
#nFeatures
#size_factors
#addPerCellQC
#addPerFeatureQC
#Mitochondrial genes ratio
#Ribosome genes ratio
#Filtering
#Cell cycle
#Cell cycle scores

#------------------------------Library---------------------------------------------------

library(SingleCellExperiment)
library(scater)
library(scran)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cowplot)
library(viridis)


#--------------------------Make directories---------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)

#dir.create("data")
dir.create("results")
dir.create("figures")
#dir.create("scripts")

#----------------------------Loading data--------------------------------------------------

#Read counts matrix and convert it to a matrix object. Change rownames to correct names
notfiltered.counts = readMM("./data/notfiltered_count_matrix.mtx")
notfiltered.rowMeatada = read.delim("./data/notfiltered_gene_row_names.txt", header = F)
rownames(notfiltered.counts) = notfiltered.rowMeatada[,2]
notfiltered.colMetadata = read.delim("./data/notfiltered_column_metadata.txt")
colnames(notfiltered.counts) = notfiltered.colMetadata[,1]

#Set SingleCellExperiment object with matrix
sce = SingleCellExperiment(assays = list(counts = notfiltered.counts)) #SingleCellExperiment object with counts matrix as assay

bioMD = read.csv("./data/snRNAseqPFC_BA10_biospecimen_metadata.csv", as.is = TRUE)

clin = read.csv("./data/snRNAseqPFC_BA10_id_mapping.csv", as.is = TRUE)[2:3]
projid = unique(clin[,1])
un_ros = unique(clin[,2])
unique_clin = data.frame(projid,un_ros)

bio_clin = merge(unique_clin,bioMD,by="projid")

supp_3 = read.csv2("./data/supp_3.csv")

bio_clin_supp = merge(bio_clin,supp_3,by="un_ros")

cell_metadata = merge(notfiltered.colMetadata,bio_clin_supp,by="projid")

#---------------------------Preparing data----------------------------------------------

#Metadata
rownames(cell_metadata) = cell_metadata[,2] #Assign cell names to row names
cell_metadata = cell_metadata[,-2] #Remove first column if cell names is the first column

colData(sce) = DataFrame(cell_metadata) #Assign metadata to colData of sce object

#-------------------------------Condition--------------------------------------------------

num.cells = length(colnames(sce)) #Number of cells in sce

#------------------------------------QC---------------------------------------------------

sce@colData$nCount = Matrix::colSums(counts(sce)) #Library size. Number of UMI counts taken across all cells
sce@colData$nFeatures = Matrix::colSums(counts(sce)>0) #Cells with any counts for each gene.
sce@colData$size_factors = scater::librarySizeFactors(sce)

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

jpeg(file="./figures/1.QC/Before_filtering/nFeatures_by_pathologygroup.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nFeatures",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_by_pathologygroup.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nCount",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_by_pathologygroup.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "percent_mito",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_vs._nFeatures.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_vs._nFeatures.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/NCells_vs_NGenes.jpeg", width=12, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=un_ros, y=log10(nCount), fill=un_ros)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by pathology group")+
  scale_fill_viridis_d(direction = -1)+
  xlab("Patients")
dev.off()

#Plot number of genes by patient
jpeg(file="./figures/1.QC/Before_filtering/NCells.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=pathology.group, fill=pathology.group)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells by pathology group")+
  scale_fill_viridis_d(direction = -1)+
  xlab("Pathology group")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/NCells_patient.jpeg", width=11, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=un_ros, fill=un_ros)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells by pathology group")+
  scale_fill_viridis_d(direction = -1)+
  xlab("Patients")
dev.off()


# Visualize the number UMIs/transcripts per cell
jpeg(file="./figures/1.QC/Before_filtering/NUMIs.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(color=pathology.group, x=nCount, fill= pathology.group))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by pathology group")+
  scale_fill_viridis_d(direction = -1)
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
abline(v=threshold)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nFeatures_AF.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "nFeatures",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCounts_AF.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "nCount",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_AF.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "percent_mito",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCounts_vs._nFeatures_AF.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,x = "nCount",y = "nFeatures",colour_by = "pathology.group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_vs_nFeatures_AF.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,x = "percent_mito",y = "nFeatures",colour_by = "pathology.group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/NCells_vs_NGenes_filtered_AF.jpeg", width=12, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=un_ros, y=nCount, fill=un_ros)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by patient after filtering")+
  xlab("Patients")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/NCells_patient.jpeg", width=11, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=un_ros, fill=un_ros)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells by pathology group")+
  scale_fill_viridis_d(direction = -1)+
  xlab("Patients")
dev.off()


#Plot number of genes by patient
jpeg(file="./figures/1.QC/After_filtering/NCells_filtered_AF.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=pathology.group, fill=pathology.group)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells by patient after filtering")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

# Visualize the number UMIs/transcripts per cell
jpeg(file="./figures/1.QC/After_filtering/NUMIs_filtered_AF.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(color=pathology.group, x=nCount, fill= pathology.group))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by patient after filtering")+
  scale_fill_viridis_d(direction = -1)
dev.off()

#-------------------------------Cell cycle--------------------------------------------------

hm.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
symb = mapIds(org.Hs.eg.db, keys=rownames(sce.filt), keytype="SYMBOL", column="ENSEMBL")
rowData(sce.filt)$SYMBOL = rownames(sce.filt)
rowData(sce.filt)$ENSEMBL = symb
assignments = cyclone(sce.filt, hm.pairs, gene.names=rowData(sce.filt)$ENSEMBL, verbose=TRUE)
sce.filt$CellCycle = assignments$phases

#Scores
sce.filt$G1.score = assignments$scores$G1
sce.filt$G2M.score = assignments$scores$G2M
sce.filt$S.score = assignments$scores$S

#Plot scores
dir.create("./figures/1.QC/Cell_cycle")

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score_vs_G1_score.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "G1.score",colour_by = "pathology.group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G1_score.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G1.score",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/S_score.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "S.score",x = "pathology.group",colour_by = "pathology.group")+
  xlab("Pathology group")+
  scale_fill_viridis_d(direction = -1)
dev.off()


#-----------------------------Save results----------------------------------------------

saveRDS(sce.filt,"results/sce_qc.rds")

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

