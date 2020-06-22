# March 2010

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Differential Expression


library(scran)
library(scater)
library(SingleCellExperiment)
library(MAST)
library(data.table)
library(edgeR)
library(dplyr)

#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/7.Differential Expression_20")

sce.hvg = readRDS("results/mark_gene_select_20.rds")


#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2 = colSums(assay(sce.hvg)>0)

jpeg(file="./figures/7.Differential Expression_20/Plot_cdr2_1.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2, y=colData(sce.hvg)$nCount)
dev.off()

colData(sce.hvg)$cngeneson = scale(cdr2)
cond = factor(colData(sce.hvg)$condition)
cond = relevel(cond, 'Control')
colData(sce.hvg)$condition=cond

sce.assay = new("SingleCellAssay",sce.hvg)
rowData(sce.assay)$primerid = rownames(sce.assay)

zlmCond = zlm(~condition + cngeneson, sce.assay,exprs_value = 'logcounts')
summaryCond = summary(zlmCond,doLRT=TRUE)
print(summaryCond, n=4)
print(summaryCond, n=4, by='D')
print(summaryCond, n=4, by='C')

summaryDt = summaryCond$datatable
chisq = summaryDt[contrast=='conditionAD' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC = summaryDt[contrast=='conditionAD' & component=='logFC',.(primerid, coef)]
fcHurdle = merge(chisq, logFC, by="primerid")

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results = data.frame(gene = fcHurdle$primerid,
                  pvalue = fcHurdle$`Pr(>Chisq)`,
                  padjusted = fcHurdle$fdr,
                  logFC = fcHurdle$coef)

#---------------------------------------Results-----------------------------------------------------------------

results = results[order(results$padjusted), ]
results_by_padj = results[order(results$padjusted), ]
results_by_absFC = results[order(abs(results$logFC),decreasing = TRUE),]

write.csv(results_by_padj,"./results/results_by_padj_20.csv")
write.csv(results_by_absFC,"./results/results_by_absFC_20.csv")

mostDE = results$gene[1:20]
results$mostDE = results$gene %in% mostDE

rownames(results) = results$gene
results$genelabels = ""
results$genelabels[1:10] =  results$gene[1:10]
results$genelabels[1:10] = rownames(results[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20/Volcano_Plot.jpeg", width=6, height=4, units="in", res=300)
ggplot(results) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE)) +
  ggtitle("Volcano plot") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------HeatMap----------------------------------------------------------------------------------

mostDE = as.character(mostDE)

library(NMF)

norm = assay(sce.hvg[mostDE, ], "logcounts")
mat = as.matrix(norm)

# heatmap
jpeg(file="./figures/7.Differential Expression_20/HeatMap.jpeg", width=6, height=4, units="in", res=300)
aheatmap(mat, annCol = colData(sce.hvg[mostDE,])$condition)
dev.off()

#---------------------------------Differential Expression (edgeR)-----------------------------------------------------

dge = DGEList(
  counts = assays(sce.hvg)$counts,
  norm.factors = rep(1, length(assays(sce.hvg)$counts[1,])),
  group = sce.hvg$condition
)

group_edgeR = factor(sce.hvg$condition)
design = model.matrix(~ group_edgeR)
dge = estimateDisp(dge, design = design, trend.method = "none")
fit = glmFit(dge, design)
lrt = glmLRT(fit)

DEedgeR <- topTags(lrt, n= Inf)
head(DEedgeR$table)
DEedgeR <- DEedgeR$table

jpeg(file="./figures/7.Differential Expression_20/edgeR.jpeg", width=6, height=4, units="in", res=300)
DEedgeR %>% ggplot(aes(x= FDR)) + geom_histogram(binwidth= 0.05) +
  ggtitle("Adjusted P Value Distribution")
dev.off()

DEedgeR %>% filter(FDR <= 0.25) %>% nrow()

deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]

jpeg(file="./figures/7.Differential Expression_20/Smear.jpeg", width=6, height=4, units="in", res=300)
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
dev.off()

#-----------------------------------------Save data----------------------------------------------------------------------

saveRDS(sce.hvg,"results/sce_diff_expr_20.rds")
saveRDS(results,"results/mostDE_list_20.rds")

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


