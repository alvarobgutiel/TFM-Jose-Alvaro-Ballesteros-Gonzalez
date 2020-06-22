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
dir.create("./figures/7.Differential Expression")

sce.hvg = readRDS("results/mark_gene_select.rds")

names(colData(sce.hvg))[names(colData(sce.hvg))=="donor_organism.development_stage.ontology_label"] = "Patients"
names(colData(sce.hvg))[names(colData(sce.hvg))=="donor_organism.sex"] = "Sex"

#---------------------------------Differential Expression (MAST)-----------------------------------------------------

cdr2 = colSums(assay(sce.hvg)>0)

jpeg(file="./figures/7.Differential Expression/Plot_cdr2.jpeg", width=6, height=4, units="in", res=300)
qplot(x=cdr2, y=colData(sce.hvg)$nCount)
dev.off()

colData(sce.hvg)$cngeneson = scale(cdr2)
cond = factor(colData(sce.hvg)$Patients)
cond = relevel(cond, 'human adult stage')
colData(sce.hvg)$Patients=cond

sce.hvg = sce.hvg[!duplicated(rownames(sce.hvg)), ]

sce.assay = new("SingleCellAssay",sce.hvg)
rowData(sce.assay)$primerid = rownames(sce.assay)

zlmCond = zlm(~ Patients + cngeneson, sce.assay,exprs_value = 'logcounts')
summaryCond = summary(zlmCond,doLRT= TRUE)
print(summaryCond, n=4)
print(summaryCond, n=4, by='D')
print(summaryCond, n=4, by='C')

summaryDt = summaryCond$datatable
chisq_1month = summaryDt[contrast=='Patients1-month-old human stage' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_1month = summaryDt[contrast=='Patients1-month-old human stage' & component=='logFC',.(primerid, coef)]
fcHurdle_1month = merge(chisq_1month, logFC_1month, by="primerid")

fcHurdle_1month[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_1month = data.frame(gene = fcHurdle_1month$primerid,
                  pvalue = fcHurdle_1month$`Pr(>Chisq)`,
                  padjusted = fcHurdle_1month$fdr,
                  logFC = fcHurdle_1month$coef)

#---------------------------------------Results-----------------------------------------------------------------

results_1month = results_1month[order(results_1month$padjusted), ]
results_by_padj_1month = results_1month[order(results_1month$padjusted), ]
results_by_absFC_1month = results_1month[order(abs(results_1month$logFC),decreasing = TRUE),]


mostDE_1month = results_1month$gene[1:20]
results_1month$mostDE_1month = results_1month$gene %in% mostDE_1month

rownames(results_1month) = results_1month$gene
results_1month$genelabels = ""
results_1month$genelabels[1:10] =  results_1month$gene[1:10]
results_1month$genelabels[1:10] = rownames(results_1month[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression/Volcano_Plot_1month.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_1month) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_1month),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_1month)) +
  ggtitle("Volcano plot - Adult vs 1 month age") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()


#-----------------------------------------------------------------------------------------------------------------------

chisq_5years = summaryDt[contrast=='Patientsyoung adult stage' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_5years = summaryDt[contrast=='Patientsyoung adult stage' & component=='logFC',.(primerid, coef)]
fcHurdle_5years = merge(chisq_5years, logFC_5years, by="primerid")

fcHurdle_5years[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_5years = data.frame(gene = fcHurdle_5years$primerid,
                            pvalue = fcHurdle_5years$`Pr(>Chisq)`,
                            padjusted = fcHurdle_5years$fdr,
                            logFC = fcHurdle_5years$coef)


results_5years = results_5years[order(results_5years$padjusted), ]
results_by_padj_5years = results_5years[order(results_5years$padjusted), ]
results_by_absFC_5years = results_5years[order(abs(results_5years$logFC),decreasing = TRUE),]


mostDE_5years = results_5years$gene[1:20]
results_5years$mostDE_5years = results_5years$gene %in% mostDE_5years

rownames(results_5years) = results_5years$gene
results_5years$genelabels = ""
results_5years$genelabels[1:10] =  results_5years$gene[1:10]
results_5years$genelabels[1:10] = rownames(results_5years[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression/Volcano_Plot_5years.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_5years) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_5years),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_5years)) +
  ggtitle("Volcano plot - Adult vs 5 years") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------

chisq_young_adult = summaryDt[contrast=='Patients5-year-old human stage' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_young_adult = summaryDt[contrast=='Patients5-year-old human stage' & component=='logFC',.(primerid, coef)]
fcHurdle_young_adult = merge(chisq_young_adult, logFC_young_adult, by="primerid")

fcHurdle_young_adult[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_young_adult = data.frame(gene = fcHurdle_young_adult$primerid,
                            pvalue = fcHurdle_young_adult$`Pr(>Chisq)`,
                            padjusted = fcHurdle_young_adult$fdr,
                            logFC = fcHurdle_young_adult$coef)


results_young_adult = results_young_adult[order(results_young_adult$padjusted), ]
results_by_padj_young_adult = results_young_adult[order(results_young_adult$padjusted), ]
results_by_absFC_young_adult = results_young_adult[order(abs(results_young_adult$logFC),decreasing = TRUE),]


mostDE_young_adult = results_young_adult$gene[1:20]
results_young_adult$mostDE_young_adult = results_young_adult$gene %in% mostDE_young_adult

rownames(results_young_adult) = results_young_adult$gene
results_young_adult$genelabels = ""
results_young_adult$genelabels[1:10] =  results_young_adult$gene[1:10]
results_young_adult$genelabels[1:10] = rownames(results_young_adult[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression/Volcano_Plot_young_adult.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_young_adult) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_young_adult),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_young_adult)) +
  ggtitle("Volcano plot - Adult vs 5 years") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------

chisq_late_adult = summaryDt[contrast=='Patientshuman late adulthood stage' & component=='H',.(primerid, `Pr(>Chisq)`)]
logFC_late_adult = summaryDt[contrast=='Patientshuman late adulthood stage' & component=='logFC',.(primerid, coef)]
fcHurdle_late_adult = merge(chisq_late_adult, logFC_late_adult, by="primerid")

fcHurdle_late_adult[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
results_late_adult = data.frame(gene = fcHurdle_late_adult$primerid,
                                 pvalue = fcHurdle_late_adult$`Pr(>Chisq)`,
                                 padjusted = fcHurdle_late_adult$fdr,
                                 logFC = fcHurdle_late_adult$coef)


results_late_adult = results_late_adult[order(results_late_adult$padjusted), ]
results_by_padj_late_adult = results_late_adult[order(results_late_adult$padjusted), ]
results_by_absFC_late_adult = results_late_adult[order(abs(results_late_adult$logFC),decreasing = TRUE),]


mostDE_late_adult = results_late_adult$gene[1:20]
results_late_adult$mostDE_late_adult = results_late_adult$gene %in% mostDE_late_adult

rownames(results_late_adult) = results_late_adult$gene
results_late_adult$genelabels = ""
results_late_adult$genelabels[1:10] =  results_late_adult$gene[1:10]
results_late_adult$genelabels[1:10] = rownames(results_late_adult[1:10,])

#---------------------------------------Volcano Plot----------------------------------------------------------------

library(ggrepel)
# volcano plot
jpeg(file="./figures/7.Differential Expression/Volcano_Plot_late_adult.jpeg", width=6, height=4, units="in", res=300)
ggplot(results_late_adult) +
  geom_text_repel(aes(x = logFC, y = -log2(pvalue), label = ifelse(genelabels != "", rownames(results_late_adult),"")),size=2) +
  geom_point(aes(x=logFC, y=-log2(pvalue), colour=mostDE_late_adult)) +
  ggtitle("Volcano plot - Adult vs 5 years") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  ylim(0, NA)
dev.off()


#-----------------------------------------Save data----------------------------------------------------------------------

saveRDS(sce.hvg,"results/sce_diff_expr.rds")
saveRDS(results_1month,"results/mostDE_list_1month.rds")
saveRDS(results_5years,"results/mostDE_list_5years.rds")
saveRDS(results_young_adult,"results/mostDE_list_young_adult.rds")
saveRDS(results_late_adult,"results/mostDE_list_late_adult.rds")

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


