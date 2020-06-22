# March 2010

# Single cell analysis - SingleCellExperiment, scater, scran
# GSE138852
# Gene Set Enrichment Analysis


library(scran)
library(scater)
library(SingleCellExperiment)
library(fgsea)
library(dplyr)
library(tibble)
library(goseq)
library(stringr)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/script")[[1]]
setwd(dir_real)
dir.create("./figures/8.Gene Set Enrichment Analysis")

sce.hvg = readRDS("results/sce_diff_expr.rds")
results = readRDS("results/mostDE_list.rds")
results_1month = readRDS("results/mostDE_list_1month.rds")
results_5years = readRDS("results/mostDE_list_5years.rds")
results_late_adult = readRDS("results/mostDE_list_late_adult.rds")
results_young_adult = readRDS("results/mostDE_list_young_adult.rds")

rownames(results_1month) = results_1month$gene
rownames(results_5years) = results_5years$gene
rownames(results_late_adult) = results_late_adult$gene
rownames(results_young_adult) = results_young_adult$gene

#---------------------------------------------fGSEA all---------------------------------------------------------------

table(results$padjusted < 0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis/logFC_barplot.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results$logFC, decreasing = T))
dev.off()

res2 = results %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks = deframe(res2)
summary(ranks)

pathways.hallmark = gmtPathways("data/msigdb.v7.1.symbols.gmt")

fgseaRes = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks, nperm=1000)
dim(fgseaRes)
colnames(fgseaRes)
table(fgseaRes$padj < 0.05)

fgseaResTidy = fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
dim(fgseaResTidy)
fgseaResTidy

jpeg(file="./figures/8.Gene Set Enrichment Analysis/Hallmark_pathways_GSEA.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy[1:30,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq all---------------------------------------------------------------

#selecciona genes con logFC diferente de 0. A continuacion ajusta los p-valores.
#y por último genera una variable de 0 y 1, donde 1 es DE y 0 NO.
genes=as.integer(p.adjust(results$pvalue[results$logFC!=0],method = "BH")<.05)
table(genes) #los resultados coinciden con la exploracion que vimos antes: 563 DE,17 no

names(genes) = rownames(results[results$logFC!=0,])

#desde goseq realizamos otro enriquecimiento:
pwf=nullp(genes,"hg19","geneSymbol",bias.data = NULL)

goResults = goseq(pwf, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults)  #6099 funciones evaluadas (solo BP)
table(goResults$ontology)
table(goResults$over_represented_pvalue <0.05 |goResults$under_represented_pvalue <0.05)
#91 funciones significativas

#representacion de las 15 funciones mas significativas:
jpeg(file="./figures/8.Gene Set Enrichment Analysis/GO_results.jpeg", width=8, height=4, units="in", res=300)
goResults %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
dev.off()

#----------------------------------------------------------------------------------------------------------------
#---------------------------------------------fGSEA 1 month---------------------------------------------------------------

table(results_1month$padjusted < 0.05)
results_ok_1month = results_1month[results_1month$padjusted < 0.05,]
results_up_1month = results_ok_1month[results_ok_1month$logFC > 0,]
results_down_1month = results_ok_1month[results_ok_1month$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis/logFC_barplot_1month.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_1month$logFC, decreasing = T))
dev.off()

res2_1month = results_1month %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_1month = deframe(res2_1month)
summary(ranks_1month)

fgseaRes_1month = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_1month, nperm=1000)
dim(fgseaRes_1month)
colnames(fgseaRes_1month)
table(fgseaRes_1month$padj < 0.05)

fgseaResTidy_1month = fgseaRes_1month %>%
  as_tibble() %>%
  arrange(desc(NES))
dim(fgseaResTidy_1month)
fgseaResTidy_1month

jpeg(file="./figures/8.Gene Set Enrichment Analysis/Hallmark_pathways_GSEA_1month1.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_1month[1:30,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA 1month") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq all---------------------------------------------------------------

genes_1month=as.integer(p.adjust(results_1month$pvalue[results_1month$logFC!=0],method = "BH")<.05)
table(genes_1month)

names(genes_1month) = rownames(results_1month[results_1month$logFC!=0,])

#desde goseq realizamos otro enriquecimiento:
pwf_1month=nullp(genes_1month,"hg19","geneSymbol",bias.data = NULL)

goResults_1month = goseq(pwf_1month, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_1month)  #6099 funciones evaluadas (solo BP)
table(goResults_1month$ontology)
table(goResults_1month$over_represented_pvalue <0.05 |goResults_1month$under_represented_pvalue <0.05)
#91 funciones significativas

#representacion de las 15 funciones mas significativas:
jpeg(file="./figures/8.Gene Set Enrichment Analysis/GO_results_1month.jpeg", width=6, height=4, units="in", res=300)
goResults_1month %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
dev.off()

#---------------------------------------------fGSEA 5 years---------------------------------------------------------------

table(results_5years$padjusted < 0.05)
results_ok_5years = results_5years[results_5years$padjusted < 0.05,]
results_up_5years = results_ok_5years[results_ok_5years$logFC > 0,]
results_down_5years = results_ok_5years[results_ok_5years$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis/logFC_barplot_5years.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_5years$logFC, decreasing = T))
dev.off()

res2_5years = results_5years %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_5years = deframe(res2_5years)
summary(ranks_5years)

fgseaRes_5years = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_5years, nperm=1000)
dim(fgseaRes_5years)
colnames(fgseaRes_5years)
table(fgseaRes_5years$padj < 0.05)

fgseaResTidy_5years = fgseaRes_5years %>%
  as_tibble() %>%
  arrange(desc(NES))
dim(fgseaResTidy_5years)
fgseaResTidy_5years

jpeg(file="./figures/8.Gene Set Enrichment Analysis/Hallmark_pathways_GSEA_5years2.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_5years[1:30,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05))+
  scale_fill_manual(values=c("#00BFC4"))+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA 5years") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq all---------------------------------------------------------------

genes_5years=as.integer(p.adjust(results_5years$pvalue[results_5years$logFC!=0],method = "BH")<.05)
table(genes_5years)

names(genes_5years) = rownames(results_5years[results_5years$logFC!=0,])

#desde goseq realizamos otro enriquecimiento:
pwf_5years=nullp(genes_5years,"hg19","geneSymbol",bias.data = NULL)

goResults_5years = goseq(pwf_5years, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_5years)  #6099 funciones evaluadas (solo BP)
table(goResults_5years$ontology)
table(goResults_5years$over_represented_pvalue <0.05 |goResults_5years$under_represented_pvalue <0.05)
#91 funciones significativas

#representacion de las 15 funciones mas significativas:
jpeg(file="./figures/8.Gene Set Enrichment Analysis/GO_results_5years.jpeg", width=6, height=4, units="in", res=300)
goResults_5years %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
dev.off()


#---------------------------------------------fGSEA young adult---------------------------------------------------------------

table(results_young_adult$padjusted < 0.05)
results_ok_young_adult = results_young_adult[results_young_adult$padjusted < 0.05,]
results_up_young_adult = results_ok_young_adult[results_ok_young_adult$logFC > 0,]
results_down_young_adult = results_ok_young_adult[results_ok_young_adult$logFC < 0,]


jpeg(file="./figures/8.Gene Set Enrichment Analysis/logFC_barplot_young_adult.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_young_adult$logFC, decreasing = T))
dev.off()

res2_young_adult = results_young_adult %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_young_adult = deframe(res2_young_adult)
summary(ranks_young_adult)

fgseaRes_young_adult = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_young_adult, nperm=1000)
dim(fgseaRes_young_adult)
colnames(fgseaRes_young_adult)
table(fgseaRes_young_adult$padj < 0.05)

fgseaResTidy_young_adult = fgseaRes_young_adult %>%
  as_tibble() %>%
  arrange(desc(NES))
dim(fgseaResTidy_young_adult)
fgseaResTidy_young_adult

jpeg(file="./figures/8.Gene Set Enrichment Analysis/Hallmark_pathways_GSEA_young_adult1.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_young_adult[1:30,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA young adult") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq all---------------------------------------------------------------

genes_young_adult=as.integer(p.adjust(results_young_adult$pvalue[results_young_adult$logFC!=0],method = "BH")<.05)
table(genes_young_adult)

names(genes_young_adult) = rownames(results_young_adult[results_young_adult$logFC!=0,])

#desde goseq realizamos otro enriquecimiento:
pwf_young_adult=nullp(genes_young_adult,"hg19","geneSymbol",bias.data = NULL)

goResults_young_adult = goseq(pwf_young_adult, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_young_adult)  #6099 funciones evaluadas (solo BP)
table(goResults_young_adult$ontology)
table(goResults_young_adult$over_represented_pvalue <0.05 |goResults_young_adult$under_represented_pvalue <0.05)
#91 funciones significativas

#representacion de las 15 funciones mas significativas:
jpeg(file="./figures/8.Gene Set Enrichment Analysis/GO_results_young_adult.jpeg", width=6, height=4, units="in", res=300)
goResults_young_adult %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
dev.off()

#---------------------------------------------fGSEA late adult---------------------------------------------------------------

table(results_late_adult$padjusted < 0.05)
results_ok_late_adult = results_late_adult[results_late_adult$padjusted < 0.05,]
results_up_late_adult = results_ok_late_adult[results_ok_late_adult$logFC > 0,]
results_down_late_adult = results_ok_late_adult[results_ok_late_adult$logFC < 0,]


jpeg(file="./figures/8.Gene Set Enrichment Analysis/logFC_barplot_late_adult.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_late_adult$logFC, decreasing = T))
dev.off()

res2_late_adult = results_late_adult %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_late_adult = deframe(res2_late_adult)
summary(ranks_late_adult)

fgseaRes_late_adult = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_late_adult, nperm=1000)
dim(fgseaRes_late_adult)
colnames(fgseaRes_late_adult)
table(fgseaRes_late_adult$padj < 0.05)

fgseaResTidy_late_adult = fgseaRes_late_adult %>%
  as_tibble() %>%
  arrange(desc(NES))
dim(fgseaResTidy_late_adult)
fgseaResTidy_late_adult

jpeg(file="./figures/8.Gene Set Enrichment Analysis/Hallmark_pathways_GSEA_late_adult1.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_late_adult[1:30,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA 1month") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq all---------------------------------------------------------------

genes_late_adult=as.integer(p.adjust(results_late_adult$pvalue[results_late_adult$logFC!=0],method = "BH")<.05)
table(genes_late_adult)

names(genes_late_adult) = rownames(results_late_adult[results_late_adult$logFC!=0,])

#desde goseq realizamos otro enriquecimiento:
pwf_late_adult=nullp(genes_late_adult,"hg19","geneSymbol",bias.data = NULL)

goResults_late_adult = goseq(pwf_late_adult, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_late_adult)  #6099 funciones evaluadas (solo BP)
table(goResults_late_adult$ontology)
table(goResults_late_adult$over_represented_pvalue <0.05 |goResults_late_adult$under_represented_pvalue <0.05)
#91 funciones significativas

#representacion de las 15 funciones mas significativas:
jpeg(file="./figures/8.Gene Set Enrichment Analysis/GO_results_late_adult.jpeg", width=6, height=4, units="in", res=300)
goResults_late_adult %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
dev.off()

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
