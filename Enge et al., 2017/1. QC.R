
# March 2010

# Single cell analysis - SingleCellExperiment
# GSE81547
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
library(RColorBrewer)


#--------------------------Make directories---------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)

#dir.create("data")
#dir.create("results")
#dir.create("figures")
#dir.create("scripts")

#----------------------------Loading data--------------------------------------------------

#Read counts matrix and convert it to a matrix object. Change rownames to correct names
counts_matrix = read.csv("./data/expression.csv")
rownames(counts_matrix) = counts_matrix[,1] #Asignar los nombres de la primera columna a las celdas
counts_matrix = counts_matrix[,-1] #Eliminamos la primera columna si contiene el nombre de las celdas
counts_matrix = as.matrix(counts_matrix) #Convertimos a matriz

genes = read.csv("./data/genes.csv")
colnames(counts_matrix) = genes$featurename[genes$featurekey == colnames(counts_matrix)]

counts_matrix = t(counts_matrix)

#Set SingleCellExperiment object with matrix
sce = SingleCellExperiment(assays = list(counts = counts_matrix)) #SingleCellExperiment object with counts matrix as assay


#---------------------------Preparing data----------------------------------------------

#Metadata
cell_metadata = read.csv("./data/cells.csv",sep=",") #Import metadata if available
rownames(cell_metadata) = cell_metadata[,1] #Assign cell names to row names
cell_metadata = cell_metadata[,-1] #Remove first column if cell names is the first column

colData(sce) = DataFrame(cell_metadata) #Assign metadata to colData of sce object

#-------------------------------Condition--------------------------------------------------

num.cells = length(colnames(sce)) #Number of cells in sce
num.cells

#Assign conditions to each cell (in this case, AD vs Control)
#for (i in 1:num.cells){
 # if (grepl("AD",sce$patient[i]) == TRUE){
  #  sce$condition[i] <- "AD"
  #}else{
   # sce$condition[i] <- "Control"
#}}

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

names(colData(sce))[names(colData(sce))=="donor_organism.development_stage.ontology_label"] = "Patients"
names(colData(sce))[names(colData(sce))=="donor_organism.sex"] = "Sex"

#Create metadata dataframe
metadata_2 = as.data.frame(sce@colData)

#Plots of QC before filtering

dir.create("./figures/1.QC/Before_filtering")

jpeg(file="./figures/1.QC/Before_filtering/density_plot_nFeatures.jpeg", width=6, height=4, units="in", res=300)
plot(density(sce$nFeatures),main="Density - total counts - Enge et al., 2017")
abline(v=200)
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nFeatures_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nFeatures",x = "Sex",colour_by = "Sex")+
  ggtitle("Detected genes per cell - By sex")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nFeatures_by_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nFeatures",x = "Patients",colour_by = "Patients")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5), legend.text = element_text(size = 8))+
  ggtitle("Detected genes per cell - By patient")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nCount",x = "Sex",colour_by = "Sex")+
  ggtitle("Total counts per cell - By sex")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_by_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nCount",x = "Patients",colour_by = "Patients")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5), legend.text = element_text(size = 8))+
  ggtitle("Total counts per cell - By patient")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "percent_mito",x = "Sex",colour_by = "Sex")+
  ggtitle("Percentage of mitochondrial genes - By sex")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_by_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "percent_mito",x = "Patients",colour_by = "Patients")+
  ggtitle("Percentage of mitochondrial genes - By patient")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5), legend.text = element_text(size = 8))
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_vs._nFeatures_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "Sex")+
  ggtitle("Number of counts vs. Number of detected genes - By sex")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/nCount_vs._nFeatures_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "Patients")+
  ggtitle("Number of counts vs. Number of detected genes - By patient")+
  theme(legend.text = element_text(size = 8))
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_vs._nFeatures_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "Sex")+
  ggtitle("Percentage of mitochondrial genes vs. Number of detected \n genes per cell - By sex")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/percent_mito_vs._nFeatures_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "Patients")+
  ggtitle("Percentage of mitochondrial genes vs. Number of detected \n genes per cell - By patient")+
  theme(legend.text = element_text(size = 8))
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/NCells_vs_NGenes_sex.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=Sex, y=log10(nCount), fill=Sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by sex")
dev.off()

#Plot number of genes by patient
jpeg(file="./figures/1.QC/Before_filtering/NCells_Sex.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=Sex, fill=Sex)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of cells by sex")
dev.off()

colfunc = colorRampPalette(c("#F7766F", "#00BBC2"))

jpeg(file="./figures/1.QC/Before_filtering/NCells_patient.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(x=Patients, fill=Patients)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of cells by patient")+
  scale_fill_manual(values = colfunc(8))
dev.off()

# Visualize the number UMIs/transcripts per cell
jpeg(file="./figures/1.QC/Before_filtering/NUMIs_sex.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(color=Sex, x=nCount, fill= Sex))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by sex")+
  theme(legend.text = element_text(size=6)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/NUMIs_patient.jpeg", width=6, height=4, units="in", res=300)
metadata_2 %>%
  ggplot(aes(color=Patients, x=nCount, fill= Patients))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by patient")+
  scale_fill_manual(values = colfunc(8))+
  theme(legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

#-------------------------------------Filtering--------------------------------------------

#Previous number of genes/cells. Data was pre-filtered.
dim(sce)

threshold = 1000 #Threshold for nFeatures (number of genes detected in each cell)

selected_c =  colnames(sce)[sce$nFeatures > threshold]
selected_f = rownames(sce)[ Matrix::rowSums(counts(sce)) > 3]
sce.filt = sce[selected_f , selected_c]

#Select Mitoratio under 0.20%
selected_mito = sce.filt$percent_mito < 0.20
# and subset the object to only keep those cells
sce.filt = sce.filt[, selected_mito]
dim(sce.filt)

#Compute the relative expression of each gene per cell
rel_expression = t(t(counts(sce.filt)) / Matrix::colSums(counts(sce.filt))) * 100
most_expressed = sort(Matrix::rowSums(rel_expression),T)[20:1] / ncol(sce.filt)

jpeg(file="./figures/1.QC/most_expressed.jpeg", width=10, height=4, units="in", res=300)
boxplot(as.matrix(t(rel_expression[names(most_expressed),])),cex=.1, las=1, xlab="% total count per cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()

#--------------------------Plot after filtering------------------------------------------------------------

dir.create("./figures/1.QC/After_filtering")

#Create metadata dataframe
metadata_3 = as.data.frame(sce.filt@colData)

jpeg(file="./figures/1.QC/After_filtering/density_plot_nFeatures.jpeg", width=6, height=4, units="in", res=300)
plot(density(sce$nFeatures),main="Density - total counts - Enge et al., 2017")
abline(v=200)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nFeatures_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nFeatures",x = "Sex",colour_by = "Sex")+
  ggtitle("Detected genes per cell - By sex")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nFeatures_by_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nFeatures",x = "Patients",colour_by = "Patients")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5), legend.text = element_text(size = 8))+
  ggtitle("Detected genes per cell - By patient")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCount_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nCount",x = "Sex",colour_by = "Sex")+
  ggtitle("Total counts per cell - By sex")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCount_by_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "nCount",x = "Patients",colour_by = "Patients")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5), legend.text = element_text(size = 8))+
  ggtitle("Total counts per cell - By patient")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_by_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "percent_mito",x = "Sex",colour_by = "Sex")+
  ggtitle("Percentage of mitochondrial genes - By sex")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_by_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,y = "percent_mito",x = "Patients",colour_by = "Patients")+
  ggtitle("Percentage of mitochondrial genes - By patient")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5), legend.text = element_text(size = 8))
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCount_vs._nFeatures_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "Sex")+
  ggtitle("Number of counts vs. Number of detected genes - By sex")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/nCount_vs._nFeatures_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "Patients")+
  ggtitle("Number of counts vs. Number of detected genes - By patient")+
  theme(legend.text = element_text(size = 8))
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_vs._nFeatures_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "Sex")+
  ggtitle("Percentage of mitochondrial genes vs. Number of detected genes per cell - By sex")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/percent_mito_vs._nFeatures_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "Patients")+
  ggtitle("Percentage of mitochondrial genes vs. Number of detected genes per cell - By patient")+
  theme(legend.text = element_text(size = 8))
dev.off()

jpeg(file="./figures/1.QC/After_filtering/NCells_vs_NGenes_sex.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=Sex, y=log10(nCount), fill=Sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes by sex")
dev.off()

#Plot number of genes by patient
jpeg(file="./figures/1.QC/After_filtering/NCells_Sex.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=Sex, fill=Sex)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of cells by sex")
dev.off()

colfunc = colorRampPalette(c("#F7766F", "#00BBC2"))

jpeg(file="./figures/1.QC/After_filtering/NCells_patient.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(x=Patients, fill=Patients)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of cells by patient")+
  scale_fill_manual(values = colfunc(8))
dev.off()

# Visualize the number UMIs/transcripts per cell
jpeg(file="./figures/1.QC/After_filtering/NUMIs_sex.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(color=Sex, x=nCount, fill= Sex))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by sex")+
  theme(legend.text = element_text(size=6)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

jpeg(file="./figures/1.QC/After_filtering/NUMIs_patient.jpeg", width=6, height=4, units="in", res=300)
metadata_3 %>%
  ggplot(aes(color=Patients, x=nCount, fill= Patients))+
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("NUMIs/transcripts per cell by patient")+
  scale_fill_manual(values = colfunc(8))+
  theme(legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
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

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score_vs_G1_score_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "G1.score",colour_by = "Sex") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("G1 score by sex")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "Sex",colour_by = "Sex") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("G2M score by sex")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G1_score_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G1.score",x = "Sex",colour_by = "Sex") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("G1 score by sex")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/S_score_sex.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "S.score",x = "Sex",colour_by = "Sex") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("S score by sex")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score_vs_G1_score_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "G1.score",colour_by = "Patients")  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("G2M score vs. G1 score by patient")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G2M_score_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G2M.score",x = "Patients",colour_by = "Patients") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("G2M score by patient")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/G1_score_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "G1.score",x = "Patients",colour_by = "Patients") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("G1 score by patient")
dev.off()

jpeg(file="./figures/1.QC/Cell_cycle/S_score_patient.jpeg", width=6, height=4, units="in", res=300)
plotColData(sce.filt,y = "S.score",x = "Patients",colour_by = "Patients") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),legend.text = element_text(size=8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("S score by patient")
dev.off()


# ------------------------------------------Library size factors ---------------------------------------------------------

lib.sce = librarySizeFactors(sce.filt)

# Examine distribution of size factors
summary(lib.sce)

jpeg(file="./figures/1.QC/log10_library_size_effect.jpeg", width=6, height=4, units="in", res=300)
hist(log10(lib.sce), xlab="Log10[Size factor]", col="grey80")
dev.off()


#-----------------------------Save results----------------------------------------------

saveRDS(sce.filt,"results/sce_qc.rds")

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

