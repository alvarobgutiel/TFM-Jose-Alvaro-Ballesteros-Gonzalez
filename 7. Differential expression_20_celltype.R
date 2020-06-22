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
cond = factor(colData(sce.hvg)$cell.type)
cond = relevel(cond, 'Astrocytes')
colData(sce.hvg)$cell.type=cond

sce.assay = new("SingleCellAssay",sce.hvg)
rowData(sce.assay)$primerid = rownames(sce.assay)

zlmCond = zlm(~cell.type + cngeneson, sce.assay,exprs_value = 'logcounts')
summaryCond = summary(zlmCond,doLRT=TRUE)
print(summaryCond, n=4)
print(summaryCond, n=4, by='D')
print(summaryCond, n=4, by='C')

summaryDt = summaryCond$datatable

chisq_endo = summaryDt[contrast=='cell.typeEndothelial Cells' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_endo = summaryDt[contrast=='cell.typeEndothelial Cells' & component=='logFC',.(primerid, coef)]
fcHurdle_endo = merge(chisq_endo, logFC_endo, by="primerid")

fcHurdle_endo[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_endo = data.frame(gene = fcHurdle_endo$primerid,
                  pvalue = fcHurdle_endo$`Pr(>Chisq)`,
                  padjusted = fcHurdle_endo$fdr,
                  logFC = fcHurdle_endo$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_endo = results_endo[order(results_endo$padjusted), ]
results_by_padj_endo = results_endo[order(results_endo$padjusted), ]
results_by_absFC_endo = results_endo[order(abs(results_endo$logFC),decreasing = TRUE),]

write.csv(results_by_padj_endo,"./results/results_by_padj_20_endo.csv")
write.csv(results_by_absFC_endo,"./results/results_by_absFC_20_endo.csv")

mostDE_endo = results_endo$gene[1:20]
results_endo$mostDE_endo = results_endo$gene %in% mostDE_endo

rownames(results_endo) = results_endo$gene
results_endo$genelabels = ""
results_endo$genelabels[1:10] =  results_endo$gene[1:10]
results_endo$genelabels[1:10] = rownames(results_endo[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20/Volcano_Plot_endo.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_endo) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_endo),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_endo)) +
  ggtitle("Volcano plot - Endothelial cells") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
#---------------------------------Differential Expression (MAST)-----------------------------------------------------

chisq_neu = summaryDt[contrast=='cell.typeNeurons' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_neu = summaryDt[contrast=='cell.typeNeurons' & component=='logFC',.(primerid, coef)]
fcHurdle_neu = merge(chisq_neu, logFC_neu, by="primerid")

fcHurdle_neu[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_neu = data.frame(gene = fcHurdle_neu$primerid,
                         pvalue = fcHurdle_neu$`Pr(>Chisq)`,
                         padjusted = fcHurdle_neu$fdr,
                         logFC = fcHurdle_neu$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_neu = results_neu[order(results_neu$padjusted), ]
results_by_padj_neu = results_neu[order(results_neu$padjusted), ]
results_by_absFC_neu = results_neu[order(abs(results_neu$logFC),decreasing = TRUE),]

write.csv(results_by_padj_neu,"./results/results_by_padj_20_neu.csv")
write.csv(results_by_absFC_neu,"./results/results_by_absFC_20_neu.csv")

mostDE_neu = results_neu$gene[1:20]
results_neu$mostDE_neu = results_neu$gene %in% mostDE_neu

rownames(results_neu) = results_neu$gene
results_neu$genelabels = ""
results_neu$genelabels[1:10] =  results_neu$gene[1:10]
results_neu$genelabels[1:10] = rownames(results_neu[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20/Volcano_Plot_neu.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_neu) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_neu),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_neu)) +
  ggtitle("Volcano plot - Neurons") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

chisq_micro = summaryDt[contrast=='cell.typeMicroglia' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_micro = summaryDt[contrast=='cell.typeMicroglia' & component=='logFC',.(primerid, coef)]
fcHurdle_micro = merge(chisq_micro, logFC_micro, by="primerid")

fcHurdle_micro[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_micro = data.frame(gene = fcHurdle_micro$primerid,
                         pvalue = fcHurdle_micro$`Pr(>Chisq)`,
                         padjusted = fcHurdle_micro$fdr,
                         logFC = fcHurdle_micro$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_micro = results_micro[order(results_micro$padjusted), ]
results_by_padj_micro = results_micro[order(results_micro$padjusted), ]
results_by_absFC_micro = results_micro[order(abs(results_micro$logFC),decreasing = TRUE),]

write.csv(results_by_padj_micro,"./results/results_by_padj_20_micro.csv")
write.csv(results_by_absFC_micro,"./results/results_by_absFC_20_micro.csv")

mostDE_micro = results_micro$gene[1:20]
results_micro$mostDE_micro = results_micro$gene %in% mostDE_micro

rownames(results_micro) = results_micro$gene
results_micro$genelabels = ""
results_micro$genelabels[1:10] =  results_micro$gene[1:10]
results_micro$genelabels[1:10] = rownames(results_micro[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20/Volcano_Plot_micro.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_micro) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_micro),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_micro)) +
  ggtitle("Volcano plot - Microglia") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

chisq_oligo = summaryDt[contrast=='cell.typeOligodendrocytes' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_oligo = summaryDt[contrast=='cell.typeOligodendrocytes' & component=='logFC',.(primerid, coef)]
fcHurdle_oligo = merge(chisq_oligo, logFC_oligo, by="primerid")

fcHurdle_oligo[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_oligo = data.frame(gene = fcHurdle_oligo$primerid,
                           pvalue = fcHurdle_oligo$`Pr(>Chisq)`,
                           padjusted = fcHurdle_oligo$fdr,
                           logFC = fcHurdle_oligo$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_oligo = results_oligo[order(results_oligo$padjusted), ]
results_by_padj_oligo = results_oligo[order(results_oligo$padjusted), ]
results_by_absFC_oligo = results_oligo[order(abs(results_oligo$logFC),decreasing = TRUE),]

write.csv(results_by_padj_oligo,"./results/results_by_padj_20_oligo.csv")
write.csv(results_by_absFC_oligo,"./results/results_by_absFC_20_oligo.csv")

mostDE_oligo = results_oligo$gene[1:20]
results_oligo$mostDE_oligo = results_oligo$gene %in% mostDE_oligo

rownames(results_oligo) = results_oligo$gene
results_oligo$genelabels = ""
results_oligo$genelabels[1:10] =  results_oligo$gene[1:10]
results_oligo$genelabels[1:10] = rownames(results_oligo[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20/Volcano_Plot_oligo.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_oligo) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_oligo),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_oligo)) +
  ggtitle("Volcano plot - Oligodendrocytes") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

chisq_oligoPC = summaryDt[contrast=='cell.typeOligodendrocyte PC' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_oligoPC = summaryDt[contrast=='cell.typeOligodendrocyte PC' & component=='logFC',.(primerid, coef)]
fcHurdle_oligoPC = merge(chisq_oligoPC, logFC_oligoPC, by="primerid")

fcHurdle_oligoPC[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_oligoPC = data.frame(gene = fcHurdle_oligoPC$primerid,
                           pvalue = fcHurdle_oligoPC$`Pr(>Chisq)`,
                           padjusted = fcHurdle_oligoPC$fdr,
                           logFC = fcHurdle_oligoPC$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_oligoPC = results_oligoPC[order(results_oligoPC$padjusted), ]
results_by_padj_oligoPC = results_oligoPC[order(results_oligoPC$padjusted), ]
results_by_absFC_oligoPC = results_oligoPC[order(abs(results_oligoPC$logFC),decreasing = TRUE),]

write.csv(results_by_padj_oligoPC,"./results/results_by_padj_20_oligoPC.csv")
write.csv(results_by_absFC_oligoPC,"./results/results_by_absFC_20_oligoPC.csv")

mostDE_oligoPC = results_oligoPC$gene[1:20]
results_oligoPC$mostDE_oligoPC = results_oligoPC$gene %in% mostDE_oligoPC

rownames(results_oligoPC) = results_oligoPC$gene
results_oligoPC$genelabels = ""
results_oligoPC$genelabels[1:10] =  results_oligoPC$gene[1:10]
results_oligoPC$genelabels[1:10] = rownames(results_oligoPC[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20/Volcano_Plot_oligoPC.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_oligoPC) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_oligoPC),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_oligoPC)) +
  ggtitle("Volcano plot - Oligodendrocytes") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()


#-----------------------------------------Save data----------------------------------------------------------------------

saveRDS(sce.hvg,"results/sce_diff_expr_20.rds")
saveRDS(results,"results/mostDE_list_20.rds")

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


