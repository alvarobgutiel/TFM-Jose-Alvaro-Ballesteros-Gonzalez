# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# Mathys et al., 2019
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
dir.create("./figures/7.Differential Expression_20_celltypes")

sce.hvg = readRDS("results/mark_gene_select_20.rds")


#---------------------------------------------Subsets---------------------------------------------------------------

sce_astro = sce.hvg[,colData(sce.hvg)$cell.type == "Astrocytes"]
sce_neu = sce.hvg[,colData(sce.hvg)$cell.type == "Neurons"]
sce_endo = sce.hvg[,colData(sce.hvg)$cell.type == "Endothelial Cells"]
sce_micro = sce.hvg[,colData(sce.hvg)$cell.type == "Microglia"]
sce_oligo = sce.hvg[,colData(sce.hvg)$cell.type == "Oligodendrocytes"]
sce_oligoPC = sce.hvg[,colData(sce.hvg)$cell.type == "Oligodendrocyte PC"]

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2_astro = colSums(assay(sce_astro)>0)

jpeg(file="./figures/7.Differential Expression_20_celltypes/Plot_cdr2_astro.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2_astro, y=colData(sce_astro)$nCount)
dev.off()

colData(sce_astro)$cngeneson = scale(cdr2_astro)
cond = factor(colData(sce_astro)$pathology.group)
cond = relevel(cond, 'no-pathology')
colData(sce_astro)$pathology.group=cond

sce.assay_astro = new("SingleCellAssay",sce_astro)
rowData(sce.assay_astro)$primerid = rownames(sce.assay_astro)

zlmCond_astro = zlm( ~ pathology.group + cngeneson, sce.assay_astro,exprs_value = 'logcounts')
summaryCond_astro = summary(zlmCond_astro,doLRT=TRUE)
print(summaryCond_astro, n=4)
print(summaryCond_astro, n=4, by='D')
print(summaryCond_astro, n=4, by='C')

summaryDt_astro = summaryCond_astro$datatable

chisq_astro_late = summaryDt_astro[contrast=='pathology.grouplate-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_astro_late = summaryDt_astro[contrast=='pathology.grouplate-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_astro_late = merge(chisq_astro_late, logFC_astro_late, by="primerid")

fcHurdle_astro_late[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_astro_late = data.frame(gene = fcHurdle_astro_late$primerid,
                  pvalue = fcHurdle_astro_late$`Pr(>Chisq)`,
                  padjusted = fcHurdle_astro_late$fdr,
                  logFC = fcHurdle_astro_late$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_astro_late = results_astro_late[order(results_astro_late$padjusted), ]
results_by_padj_astro_late = results_astro_late[order(results_astro_late$padjusted), ]
results_by_absFC_astro_late = results_astro_late[order(abs(results_astro_late$logFC),decreasing = TRUE),]

write.csv(results_by_padj_astro_late,"./results/results_by_padj_20_astroAD_late.csv")
write.csv(results_by_absFC_astro_late,"./results/results_by_absFC_20_astroAD_late.csv")

mostDE_astro_late = results_astro_late$gene[1:20]
results_astro_late$mostDE_astro_late = results_astro_late$gene %in% mostDE_astro_late

rownames(results_astro_late) = results_astro_late$gene
results_astro_late$genelabels = ""
results_astro_late$genelabels[1:10] =  results_astro_late$gene[1:10]
results_astro_late$genelabels[1:10] = rownames(results_astro_late[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_astro_late.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_astro_late) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_astro_late),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_astro_late)) +
  ggtitle("Volcano plot - AD Astrocytes Late AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#------------------------------------------------------------------------------------------------

chisq_astro_early = summaryDt_astro[contrast=='pathology.groupearly-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_astro_early = summaryDt_astro[contrast=='pathology.groupearly-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_astro_early = merge(chisq_astro_early, logFC_astro_early, by="primerid")

fcHurdle_astro_early[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_astro_early = data.frame(gene = fcHurdle_astro_early$primerid,
                                pvalue = fcHurdle_astro_early$`Pr(>Chisq)`,
                                padjusted = fcHurdle_astro_early$fdr,
                                logFC = fcHurdle_astro_early$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_astro_early = results_astro_early[order(results_astro_early$padjusted), ]
results_by_padj_astro_early = results_astro_early[order(results_astro_early$padjusted), ]
results_by_absFC_astro_early = results_astro_early[order(abs(results_astro_early$logFC),decreasing = TRUE),]

write.csv(results_by_padj_astro_early,"./results/results_by_padj_20_astroAD_late.csv")
write.csv(results_by_absFC_astro_early,"./results/results_by_absFC_20_astroAD_late.csv")

mostDE_astro_early = results_astro_early$gene[1:20]
results_astro_early$mostDE_astro_early = results_astro_early$gene %in% mostDE_astro_early

rownames(results_astro_early) = results_astro_early$gene
results_astro_early$genelabels = ""
results_astro_early$genelabels[1:10] =  results_astro_early$gene[1:10]
results_astro_early$genelabels[1:10] = rownames(results_astro_early[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_astro_early.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_astro_early) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_astro_early),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_astro_early)) +
  ggtitle("Volcano plot - AD Astrocytes Early AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2_neu = colSums(assay(sce_neu)>0)

jpeg(file="./figures/7.Differential Expression_20_celltypes/Plot_cdr2_neu.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2_neu, y=colData(sce_neu)$nCount)
dev.off()

colData(sce_neu)$cngeneson = scale(cdr2_neu)
cond = factor(colData(sce_neu)$pathology.group)
cond = relevel(cond, 'no-pathology')
colData(sce_neu)$pathology.group=cond

sce.assay_neu = new("SingleCellAssay",sce_neu)
rowData(sce.assay_neu)$primerid = rownames(sce.assay_neu)

zlmCond_neu = zlm( ~ pathology.group + cngeneson, sce.assay_neu,exprs_value = 'logcounts')
summaryCond_neu = summary(zlmCond_neu,doLRT=TRUE)
print(summaryCond_neu, n=4)
print(summaryCond_neu, n=4, by='D')
print(summaryCond_neu, n=4, by='C')

summaryDt_neu = summaryCond_neu$datatable

chisq_neu_late = summaryDt_neu[contrast=='pathology.grouplate-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_neu_late = summaryDt_neu[contrast=='pathology.grouplate-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_neu_late = merge(chisq_neu_late, logFC_neu_late, by="primerid")

fcHurdle_neu_late[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_neu_late = data.frame(gene = fcHurdle_neu_late$primerid,
                                pvalue = fcHurdle_neu_late$`Pr(>Chisq)`,
                                padjusted = fcHurdle_neu_late$fdr,
                                logFC = fcHurdle_neu_late$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_neu_late = results_neu_late[order(results_neu_late$padjusted), ]
results_by_padj_neu_late = results_neu_late[order(results_neu_late$padjusted), ]
results_by_absFC_neu_late = results_neu_late[order(abs(results_neu_late$logFC),decreasing = TRUE),]

write.csv(results_by_padj_neu_late,"./results/results_by_padj_20_neuAD_late.csv")
write.csv(results_by_absFC_neu_late,"./results/results_by_absFC_20_neuAD_late.csv")

mostDE_neu_late = results_neu_late$gene[1:20]
results_neu_late$mostDE_neu_late = results_neu_late$gene %in% mostDE_neu_late

rownames(results_neu_late) = results_neu_late$gene
results_neu_late$genelabels = ""
results_neu_late$genelabels[1:10] =  results_neu_late$gene[1:10]
results_neu_late$genelabels[1:10] = rownames(results_neu_late[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_neu_late.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_neu_late) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_neu_late),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_neu_late)) +
  ggtitle("Volcano plot - AD Neurons Late AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#------------------------------------------------------------------------------------------------

chisq_neu_early = summaryDt_neu[contrast=='pathology.groupearly-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_neu_early = summaryDt_neu[contrast=='pathology.groupearly-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_neu_early = merge(chisq_neu_early, logFC_neu_early, by="primerid")

fcHurdle_neu_early[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_neu_early = data.frame(gene = fcHurdle_neu_early$primerid,
                                 pvalue = fcHurdle_neu_early$`Pr(>Chisq)`,
                                 padjusted = fcHurdle_neu_early$fdr,
                                 logFC = fcHurdle_neu_early$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_neu_early = results_neu_early[order(results_neu_early$padjusted), ]
results_by_padj_neu_early = results_neu_early[order(results_neu_early$padjusted), ]
results_by_absFC_neu_early = results_neu_early[order(abs(results_neu_early$logFC),decreasing = TRUE),]

write.csv(results_by_padj_neu_early,"./results/results_by_padj_20_neuAD_late.csv")
write.csv(results_by_absFC_neu_early,"./results/results_by_absFC_20_neuAD_late.csv")

mostDE_neu_early = results_neu_early$gene[1:20]
results_neu_early$mostDE_neu_early = results_neu_early$gene %in% mostDE_neu_early

rownames(results_neu_early) = results_neu_early$gene
results_neu_early$genelabels = ""
results_neu_early$genelabels[1:10] =  results_neu_early$gene[1:10]
results_neu_early$genelabels[1:10] = rownames(results_neu_early[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_neu_early.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_neu_early) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_neu_early),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_neu_early)) +
  ggtitle("Volcano plot - AD Neurons Early AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2_endo = colSums(assay(sce_endo)>0)

jpeg(file="./figures/7.Differential Expression_20_celltypes/Plot_cdr2_endo.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2_endo, y=colData(sce_endo)$nCount)
dev.off()

colData(sce_endo)$cngeneson = scale(cdr2_endo)
cond = factor(colData(sce_endo)$pathology.group)
cond = relevel(cond, 'no-pathology')
colData(sce_endo)$pathology.group=cond

sce.assay_endo = new("SingleCellAssay",sce_endo)
rowData(sce.assay_endo)$primerid = rownames(sce.assay_endo)

zlmCond_endo = zlm( ~ pathology.group + cngeneson, sce.assay_endo,exprs_value = 'logcounts')
summaryCond_endo = summary(zlmCond_endo,doLRT=TRUE)
print(summaryCond_endo, n=4)
print(summaryCond_endo, n=4, by='D')
print(summaryCond_endo, n=4, by='C')

summaryDt_endo = summaryCond_endo$datatable

chisq_endo_late = summaryDt_endo[contrast=='pathology.grouplate-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_endo_late = summaryDt_endo[contrast=='pathology.grouplate-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_endo_late = merge(chisq_endo_late, logFC_endo_late, by="primerid")

fcHurdle_endo_late[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_endo_late = data.frame(gene = fcHurdle_endo_late$primerid,
                              pvalue = fcHurdle_endo_late$`Pr(>Chisq)`,
                              padjusted = fcHurdle_endo_late$fdr,
                              logFC = fcHurdle_endo_late$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_endo_late = results_endo_late[order(results_endo_late$padjusted), ]
results_by_padj_endo_late = results_endo_late[order(results_endo_late$padjusted), ]
results_by_absFC_endo_late = results_endo_late[order(abs(results_endo_late$logFC),decreasing = TRUE),]

write.csv(results_by_padj_endo_late,"./results/results_by_padj_20_endoAD_late.csv")
write.csv(results_by_absFC_endo_late,"./results/results_by_absFC_20_endoAD_late.csv")

mostDE_endo_late = results_endo_late$gene[1:20]
results_endo_late$mostDE_endo_late = results_endo_late$gene %in% mostDE_endo_late

rownames(results_endo_late) = results_endo_late$gene
results_endo_late$genelabels = ""
results_endo_late$genelabels[1:10] =  results_endo_late$gene[1:10]
results_endo_late$genelabels[1:10] = rownames(results_endo_late[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_endo_late.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_endo_late) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_endo_late),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_endo_late)) +
  ggtitle("Volcano plot - AD Endothelial Cells Late AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#------------------------------------------------------------------------------------------------

chisq_endo_early = summaryDt_endo[contrast=='pathology.groupearly-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_endo_early = summaryDt_endo[contrast=='pathology.groupearly-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_endo_early = merge(chisq_endo_early, logFC_endo_early, by="primerid")

fcHurdle_endo_early[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_endo_early = data.frame(gene = fcHurdle_endo_early$primerid,
                               pvalue = fcHurdle_endo_early$`Pr(>Chisq)`,
                               padjusted = fcHurdle_endo_early$fdr,
                               logFC = fcHurdle_endo_early$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_endo_early = results_endo_early[order(results_endo_early$padjusted), ]
results_by_padj_endo_early = results_endo_early[order(results_endo_early$padjusted), ]
results_by_absFC_endo_early = results_endo_early[order(abs(results_endo_early$logFC),decreasing = TRUE),]

write.csv(results_by_padj_endo_early,"./results/results_by_padj_20_endoAD_late.csv")
write.csv(results_by_absFC_endo_early,"./results/results_by_absFC_20_endoAD_late.csv")

mostDE_endo_early = results_endo_early$gene[1:20]
results_endo_early$mostDE_endo_early = results_endo_early$gene %in% mostDE_endo_early

rownames(results_endo_early) = results_endo_early$gene
results_endo_early$genelabels = ""
results_endo_early$genelabels[1:10] =  results_endo_early$gene[1:10]
results_endo_early$genelabels[1:10] = rownames(results_endo_early[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_endo_early.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_endo_early) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_endo_early),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_endo_early)) +
  ggtitle("Volcano plot - AD Endothelial Cells Early AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2_micro = colSums(assay(sce_micro)>0)

jpeg(file="./figures/7.Differential Expression_20_celltypes/Plot_cdr2_micro.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2_micro, y=colData(sce_micro)$nCount)
dev.off()

colData(sce_micro)$cngeneson = scale(cdr2_micro)
cond = factor(colData(sce_micro)$pathology.group)
cond = relevel(cond, 'no-pathology')
colData(sce_micro)$pathology.group=cond

sce.assay_micro = new("SingleCellAssay",sce_micro)
rowData(sce.assay_micro)$primerid = rownames(sce.assay_micro)

zlmCond_micro = zlm( ~ pathology.group + cngeneson, sce.assay_micro,exprs_value = 'logcounts')
summaryCond_micro = summary(zlmCond_micro,doLRT=TRUE)
print(summaryCond_micro, n=4)
print(summaryCond_micro, n=4, by='D')
print(summaryCond_micro, n=4, by='C')

summaryDt_micro = summaryCond_micro$datatable

chisq_micro_late = summaryDt_micro[contrast=='pathology.grouplate-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_micro_late = summaryDt_micro[contrast=='pathology.grouplate-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_micro_late = merge(chisq_micro_late, logFC_micro_late, by="primerid")

fcHurdle_micro_late[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_micro_late = data.frame(gene = fcHurdle_micro_late$primerid,
                               pvalue = fcHurdle_micro_late$`Pr(>Chisq)`,
                               padjusted = fcHurdle_micro_late$fdr,
                               logFC = fcHurdle_micro_late$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_micro_late = results_micro_late[order(results_micro_late$padjusted), ]
results_by_padj_micro_late = results_micro_late[order(results_micro_late$padjusted), ]
results_by_absFC_micro_late = results_micro_late[order(abs(results_micro_late$logFC),decreasing = TRUE),]

write.csv(results_by_padj_micro_late,"./results/results_by_padj_20_microAD_late.csv")
write.csv(results_by_absFC_micro_late,"./results/results_by_absFC_20_microAD_late.csv")

mostDE_micro_late = results_micro_late$gene[1:20]
results_micro_late$mostDE_micro_late = results_micro_late$gene %in% mostDE_micro_late

rownames(results_micro_late) = results_micro_late$gene
results_micro_late$genelabels = ""
results_micro_late$genelabels[1:10] =  results_micro_late$gene[1:10]
results_micro_late$genelabels[1:10] = rownames(results_micro_late[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_micro_late.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_micro_late) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_micro_late),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_micro_late)) +
  ggtitle("Volcano plot - AD Microglia Late AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#------------------------------------------------------------------------------------------------

chisq_micro_early = summaryDt_micro[contrast=='pathology.groupearly-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_micro_early = summaryDt_micro[contrast=='pathology.groupearly-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_micro_early = merge(chisq_micro_early, logFC_micro_early, by="primerid")

fcHurdle_micro_early[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_micro_early = data.frame(gene = fcHurdle_micro_early$primerid,
                                pvalue = fcHurdle_micro_early$`Pr(>Chisq)`,
                                padjusted = fcHurdle_micro_early$fdr,
                                logFC = fcHurdle_micro_early$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_micro_early = results_micro_early[order(results_micro_early$padjusted), ]
results_by_padj_micro_early = results_micro_early[order(results_micro_early$padjusted), ]
results_by_absFC_micro_early = results_micro_early[order(abs(results_micro_early$logFC),decreasing = TRUE),]

write.csv(results_by_padj_micro_early,"./results/results_by_padj_20_microAD_late.csv")
write.csv(results_by_absFC_micro_early,"./results/results_by_absFC_20_microAD_late.csv")

mostDE_micro_early = results_micro_early$gene[1:20]
results_micro_early$mostDE_micro_early = results_micro_early$gene %in% mostDE_micro_early

rownames(results_micro_early) = results_micro_early$gene
results_micro_early$genelabels = ""
results_micro_early$genelabels[1:10] =  results_micro_early$gene[1:10]
results_micro_early$genelabels[1:10] = rownames(results_micro_early[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_micro_early.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_micro_early) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_micro_early),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_micro_early)) +
  ggtitle("Volcano plot - AD Microglia Early AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2_oligo = colSums(assay(sce_oligo)>0)

jpeg(file="./figures/7.Differential Expression_20_celltypes/Plot_cdr2_oligo.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2_oligo, y=colData(sce_oligo)$nCount)
dev.off()

colData(sce_oligo)$cngeneson = scale(cdr2_oligo)
cond = factor(colData(sce_oligo)$pathology.group)
cond = relevel(cond, 'no-pathology')
colData(sce_oligo)$pathology.group=cond

sce.assay_oligo = new("SingleCellAssay",sce_oligo)
rowData(sce.assay_oligo)$primerid = rownames(sce.assay_oligo)

zlmCond_oligo = zlm( ~ pathology.group + cngeneson, sce.assay_oligo,exprs_value = 'logcounts')
summaryCond_oligo = summary(zlmCond_oligo,doLRT=TRUE)
print(summaryCond_oligo, n=4)
print(summaryCond_oligo, n=4, by='D')
print(summaryCond_oligo, n=4, by='C')

summaryDt_oligo = summaryCond_oligo$datatable

chisq_oligo_late = summaryDt_oligo[contrast=='pathology.grouplate-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_oligo_late = summaryDt_oligo[contrast=='pathology.grouplate-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_oligo_late = merge(chisq_oligo_late, logFC_oligo_late, by="primerid")

fcHurdle_oligo_late[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_oligo_late = data.frame(gene = fcHurdle_oligo_late$primerid,
                                pvalue = fcHurdle_oligo_late$`Pr(>Chisq)`,
                                padjusted = fcHurdle_oligo_late$fdr,
                                logFC = fcHurdle_oligo_late$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_oligo_late = results_oligo_late[order(results_oligo_late$padjusted), ]
results_by_padj_oligo_late = results_oligo_late[order(results_oligo_late$padjusted), ]
results_by_absFC_oligo_late = results_oligo_late[order(abs(results_oligo_late$logFC),decreasing = TRUE),]

write.csv(results_by_padj_oligo_late,"./results/results_by_padj_20_oligoAD_late.csv")
write.csv(results_by_absFC_oligo_late,"./results/results_by_absFC_20_oligoAD_late.csv")

mostDE_oligo_late = results_oligo_late$gene[1:20]
results_oligo_late$mostDE_oligo_late = results_oligo_late$gene %in% mostDE_oligo_late

rownames(results_oligo_late) = results_oligo_late$gene
results_oligo_late$genelabels = ""
results_oligo_late$genelabels[1:10] =  results_oligo_late$gene[1:10]
results_oligo_late$genelabels[1:10] = rownames(results_oligo_late[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_oligo_late.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_oligo_late) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_oligo_late),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_oligo_late)) +
  ggtitle("Volcano plot - AD Oligodencrocytes Late AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#------------------------------------------------------------------------------------------------

chisq_oligo_early = summaryDt_oligo[contrast=='pathology.groupearly-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_oligo_early = summaryDt_oligo[contrast=='pathology.groupearly-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_oligo_early = merge(chisq_oligo_early, logFC_oligo_early, by="primerid")

fcHurdle_oligo_early[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_oligo_early = data.frame(gene = fcHurdle_oligo_early$primerid,
                                 pvalue = fcHurdle_oligo_early$`Pr(>Chisq)`,
                                 padjusted = fcHurdle_oligo_early$fdr,
                                 logFC = fcHurdle_oligo_early$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_oligo_early = results_oligo_early[order(results_oligo_early$padjusted), ]
results_by_padj_oligo_early = results_oligo_early[order(results_oligo_early$padjusted), ]
results_by_absFC_oligo_early = results_oligo_early[order(abs(results_oligo_early$logFC),decreasing = TRUE),]

write.csv(results_by_padj_oligo_early,"./results/results_by_padj_20_oligoAD_late.csv")
write.csv(results_by_absFC_oligo_early,"./results/results_by_absFC_20_oligoAD_late.csv")

mostDE_oligo_early = results_oligo_early$gene[1:20]
results_oligo_early$mostDE_oligo_early = results_oligo_early$gene %in% mostDE_oligo_early

rownames(results_oligo_early) = results_oligo_early$gene
results_oligo_early$genelabels = ""
results_oligo_early$genelabels[1:10] =  results_oligo_early$gene[1:10]
results_oligo_early$genelabels[1:10] = rownames(results_oligo_early[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_oligo_early.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_oligo_early) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_oligo_early),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_oligo_early)) +
  ggtitle("Volcano plot - AD Oligodencrocytes Early AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

cdr2_oligoPC = colSums(assay(sce_oligoPC)>0)

jpeg(file="./figures/7.Differential Expression_20_celltypes/Plot_cdr2_oligoPC.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2_oligoPC, y=colData(sce_oligoPC)$nCount)
dev.off()

colData(sce_oligoPC)$cngeneson = scale(cdr2_oligoPC)
cond = factor(colData(sce_oligoPC)$pathology.group)
cond = relevel(cond, 'no-pathology')
colData(sce_oligoPC)$pathology.group=cond

sce.assay_oligoPC = new("SingleCellAssay",sce_oligoPC)
rowData(sce.assay_oligoPC)$primerid = rownames(sce.assay_oligoPC)

zlmCond_oligoPC = zlm( ~ pathology.group + cngeneson, sce.assay_oligoPC,exprs_value = 'logcounts')
summaryCond_oligoPC = summary(zlmCond_oligoPC,doLRT=TRUE)
print(summaryCond_oligoPC, n=4)
print(summaryCond_oligoPC, n=4, by='D')
print(summaryCond_oligoPC, n=4, by='C')

summaryDt_oligoPC = summaryCond_oligoPC$datatable

chisq_oligoPC_late = summaryDt_oligoPC[contrast=='pathology.grouplate-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_oligoPC_late = summaryDt_oligoPC[contrast=='pathology.grouplate-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_oligoPC_late = merge(chisq_oligoPC_late, logFC_oligoPC_late, by="primerid")

fcHurdle_oligoPC_late[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_oligoPC_late = data.frame(gene = fcHurdle_oligoPC_late$primerid,
                                pvalue = fcHurdle_oligoPC_late$`Pr(>Chisq)`,
                                padjusted = fcHurdle_oligoPC_late$fdr,
                                logFC = fcHurdle_oligoPC_late$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_oligoPC_late = results_oligoPC_late[order(results_oligoPC_late$padjusted), ]
results_by_padj_oligoPC_late = results_oligoPC_late[order(results_oligoPC_late$padjusted), ]
results_by_absFC_oligoPC_late = results_oligoPC_late[order(abs(results_oligoPC_late$logFC),decreasing = TRUE),]

write.csv(results_by_padj_oligoPC_late,"./results/results_by_padj_20_oligoPCAD_late.csv")
write.csv(results_by_absFC_oligoPC_late,"./results/results_by_absFC_20_oligoPCAD_late.csv")

mostDE_oligoPC_late = results_oligoPC_late$gene[1:20]
results_oligoPC_late$mostDE_oligoPC_late = results_oligoPC_late$gene %in% mostDE_oligoPC_late

rownames(results_oligoPC_late) = results_oligoPC_late$gene
results_oligoPC_late$genelabels = ""
results_oligoPC_late$genelabels[1:10] =  results_oligoPC_late$gene[1:10]
results_oligoPC_late$genelabels[1:10] = rownames(results_oligoPC_late[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_oligoPC_late.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_oligoPC_late) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_oligoPC_late),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_oligoPC_late)) +
  ggtitle("Volcano plot - AD Oligodencrocytes Late AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#------------------------------------------------------------------------------------------------

chisq_oligoPC_early = summaryDt_oligoPC[contrast=='pathology.groupearly-pathology' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_oligoPC_early = summaryDt_oligoPC[contrast=='pathology.groupearly-pathology' & component=='logFC',.(primerid, coef)]
fcHurdle_oligoPC_early = merge(chisq_oligoPC_early, logFC_oligoPC_early, by="primerid")

fcHurdle_oligoPC_early[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_oligoPC_early = data.frame(gene = fcHurdle_oligoPC_early$primerid,
                                 pvalue = fcHurdle_oligoPC_early$`Pr(>Chisq)`,
                                 padjusted = fcHurdle_oligoPC_early$fdr,
                                 logFC = fcHurdle_oligoPC_early$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_oligoPC_early = results_oligoPC_early[order(results_oligoPC_early$padjusted), ]
results_by_padj_oligoPC_early = results_oligoPC_early[order(results_oligoPC_early$padjusted), ]
results_by_absFC_oligoPC_early = results_oligoPC_early[order(abs(results_oligoPC_early$logFC),decreasing = TRUE),]

write.csv(results_by_padj_oligoPC_early,"./results/results_by_padj_20_oligoPCAD_late.csv")
write.csv(results_by_absFC_oligoPC_early,"./results/results_by_absFC_20_oligoPCAD_late.csv")

mostDE_oligoPC_early = results_oligoPC_early$gene[1:20]
results_oligoPC_early$mostDE_oligoPC_early = results_oligoPC_early$gene %in% mostDE_oligoPC_early

rownames(results_oligoPC_early) = results_oligoPC_early$gene
results_oligoPC_early$genelabels = ""
results_oligoPC_early$genelabels[1:10] =  results_oligoPC_early$gene[1:10]
results_oligoPC_early$genelabels[1:10] = rownames(results_oligoPC_early[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression_20_celltypes/Volcano_Plot_oligoPC_early.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_oligoPC_early) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_oligoPC_early),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_oligoPC_early)) +
  ggtitle("Volcano plot - AD Oligodencrocytes Early AD") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


#-----------------------------------------Save data----------------------------------------------------------------------

saveRDS(sce.hvg,"results/sce_diff_expr_20.rds")
saveRDS(results_astro_early,"results/mostDE_list_20_astro_early.rds")
saveRDS(results_micro_early,"results/mostDE_list_20_micro_early.rds")
saveRDS(results_endo_early,"results/mostDE_list_20_endo_early.rds")
saveRDS(results_neu_early,"results/mostDE_list_20_neu_early.rds")
saveRDS(results_oligo_early,"results/mostDE_list_20_oligo_early.rds")
saveRDS(results_oligoPC_early,"results/mostDE_list_20_oligoPC_early.rds")

saveRDS(results_astro_late,"results/mostDE_list_20_astro_late.rds")
saveRDS(results_micro_late,"results/mostDE_list_20_micro_late.rds")
saveRDS(results_endo_late,"results/mostDE_list_20_endo_late.rds")
saveRDS(results_neu_late,"results/mostDE_list_20_neu_late.rds")
saveRDS(results_oligo_late,"results/mostDE_list_20_oligo_late.rds")
saveRDS(results_oligoPC_late,"results/mostDE_list_20_oligoPC_late.rds")

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


