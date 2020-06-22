# March 2020

# Single cell analysis - SingleCellExperiment, scater, scran
# Mathys et al., 2019
# Gene Set Enrichment Analysis


library(scran)
library(scater)
library(SingleCellExperiment)
library(fgsea)
library(dplyr)
library(tibble)
library(goseq)
library(stringr)
library(EnsDb.Hsapiens.v79)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/8.Gene Set Enrichment Analysis_20_cell_type_early")

sce.hvg = readRDS("results/sce_diff_expr.rds")
results_astro_early = readRDS("results/mostDE_list_20_astro_early.rds")
results_micro_early = readRDS("results/mostDE_list_20_micro_early.rds")
results_endo_early = readRDS("results/mostDE_list_20_endo_early.rds")
results_neu_early = readRDS("results/mostDE_list_20_neu_early.rds")
results_oligo_early = readRDS("results/mostDE_list_20_oligo_early.rds")
results_oligoPC_early = readRDS("results/mostDE_list_20_oligoPC_early.rds")

#---------------------------------------------fGSEA astro early---------------------------------------------------------------

table(results_astro_early$padjusted < 0.05)
results_ok_astro_early = results_astro_early[results_astro_early$padjusted < 0.05,]
results_up_astro_early = results_ok_astro_early[results_ok_astro_early$logFC > 0,]
results_down_astro_early = results_ok_astro_early[results_ok_astro_early$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/logFC_barplot_results_astro_early.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_astro_early$logFC, decreasing = T))
dev.off()

res2_astro = results_astro_early %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene)

ranks_astro = deframe(res2_astro)
summary(ranks_astro)

pathways.hallmark = gmtPathways("data/msigdb.v7.1.symbols.gmt")

fgseaRes_astro = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_astro)
dim(fgseaRes_astro)
colnames(fgseaRes_astro)
table(fgseaRes_astro$padj < 0.05)

fgseaResTidy_astro = fgseaRes_astro %>%
  as_tibble() %>%
  arrange(padj)

dim(fgseaResTidy_astro)
fgseaResTidy_astro


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/Hallmark_pathways_GSEA_astro.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_astro[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq astro early---------------------------------------------------------------


genes_astro=as.integer(p.adjust(results_astro_early$pvalue[results_astro_early$logFC!=0],method = "BH")<.05)

names(genes_astro) = rownames(results_astro_early[results_astro_early$logFC!=0,])
table(genes_astro)

pwf_astro=nullp(genes_astro,"hg19","geneSymbol",bias.data = NULL)

goResults_astro = goseq(pwf_astro, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_astro)  #6099 funciones evaluadas (solo BP)
table(goResults_astro$ontology)
table(goResults_astro$over_represented_pvalue <0.05 |goResults_astro$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/GO_results_astro.jpeg", width=8, height=6, units="in", res=300)
goResults_astro %>% 
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

#---------------------------------------------fGSEA micro early---------------------------------------------------------------

table(results_micro_early$padjusted < 0.05)
results_ok_micro_early = results_micro_early[results_micro_early$padjusted < 0.05,]
results_up_micro_early = results_ok_micro_early[results_ok_micro_early$logFC > 0,]
results_down_micro_early = results_ok_micro_early[results_ok_micro_early$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/logFC_barplot_results_micro_early.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_micro_early$logFC, decreasing = T))
dev.off()

res2_micro = results_micro_early %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene)

ranks_micro = deframe(res2_micro)
summary(ranks_micro)

fgseaRes_micro = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_micro)
dim(fgseaRes_micro)
colnames(fgseaRes_micro)
table(fgseaRes_micro$padj < 0.05)

fgseaResTidy_micro = fgseaRes_micro %>%
  as_tibble() %>%
  arrange(padj)

dim(fgseaResTidy_micro)
fgseaResTidy_micro


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/Hallmark_pathways_GSEA_micro.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_micro[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +  scale_fill_manual(values=c("#00BFC4"))+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq micro early---------------------------------------------------------------


genes_micro=as.integer(p.adjust(results_micro_early$pvalue[results_micro_early$logFC!=0],method = "BH")<.05)

names(genes_micro) = rownames(results_micro_early[results_micro_early$logFC!=0,])
table(genes_micro)

pwf_micro=nullp(genes_micro,"hg19","geneSymbol",bias.data = NULL)

goResults_micro = goseq(pwf_micro, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_micro)  #6099 funciones evaluadas (solo BP)
table(goResults_micro$ontology)
table(goResults_micro$over_represented_pvalue <0.05 |goResults_micro$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/GO_results_micro.jpeg", width=6, height=4, units="in", res=300)
goResults_micro %>% 
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

#---------------------------------------------fGSEA endo early---------------------------------------------------------------

table(results_endo_early$padjusted < 0.05)
results_ok_endo_early = results_endo_early[results_endo_early$padjusted < 0.05,]
results_up_endo_early = results_ok_endo_early[results_ok_endo_early$logFC > 0,]
results_down_endo_early = results_ok_endo_early[results_ok_endo_early$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/logFC_barplot_results_endo_early.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_endo_early$logFC, decreasing = T))
dev.off()

res2_endo = results_endo_early %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene)

ranks_endo = deframe(res2_endo)
summary(ranks_endo)

fgseaRes_endo = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_endo)
dim(fgseaRes_endo)
colnames(fgseaRes_endo)
table(fgseaRes_endo$padj < 0.05)

fgseaResTidy_endo = fgseaRes_endo %>%
  as_tibble() %>%
  arrange(padj)

dim(fgseaResTidy_endo)
fgseaResTidy_endo


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/Hallmark_pathways_GSEA_endo.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_endo[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq endo early---------------------------------------------------------------


genes_endo=as.integer(p.adjust(results_endo_early$pvalue[results_endo_early$logFC!=0],method = "BH")<.05)

names(genes_endo) = rownames(results_endo_early[results_endo_early$logFC!=0,])
table(genes_endo)

pwf_endo=nullp(genes_endo,"hg19","geneSymbol",bias.data = NULL)

goResults_endo = goseq(pwf_endo, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_endo)  #6099 funciones evaluadas (solo BP)
table(goResults_endo$ontology)
table(goResults_endo$over_represented_pvalue <0.05 |goResults_endo$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/GO_results_endo.jpeg", width=6, height=4, units="in", res=300)
goResults_endo %>% 
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

#---------------------------------------------fGSEA neu early---------------------------------------------------------------

table(results_neu_early$padjusted < 0.05)
results_ok_neu_early = results_neu_early[results_neu_early$padjusted < 0.05,]
results_up_neu_early = results_ok_neu_early[results_ok_neu_early$logFC > 0,]
results_down_neu_early = results_ok_neu_early[results_ok_neu_early$logFC < 0,]


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/logFC_barplot_results_neu_early.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_neu_early$logFC, decreasing = T))
dev.off()

res2_neu = results_neu_early %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene)

ranks_neu = deframe(res2_neu)
summary(ranks_neu)

fgseaRes_neu = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_neu)
dim(fgseaRes_neu)
colnames(fgseaRes_neu)
table(fgseaRes_neu$padj < 0.05)

fgseaResTidy_neu = fgseaRes_neu %>%
  as_tibble() %>%
  arrange(padj)

dim(fgseaResTidy_neu)
fgseaResTidy_neu


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/Hallmark_pathways_GSEA_neu.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_neu[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +  scale_fill_manual(values=c("#00BFC4"))+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq neu early---------------------------------------------------------------

genes_neu=as.integer(p.adjust(results_neu_early$pvalue[results_neu_early$logFC!=0],method = "BH")<.05)

names(genes_neu) = rownames(results_neu_early[results_neu_early$logFC!=0,])
table(genes_neu)

pwf_neu=nullp(genes_neu,"hg19","geneSymbol",bias.data = NULL)

goResults_neu = goseq(pwf_neu, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_neu)  #6099 funciones evaluadas (solo BP)
table(goResults_neu$ontology)
table(goResults_neu$over_represented_pvalue <0.05 |goResults_neu$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/GO_results_neu.jpeg", width=6, height=4, units="in", res=300)
goResults_neu %>% 
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

#---------------------------------------------fGSEA oligo early---------------------------------------------------------------

table(results_oligo_early$padjusted < 0.05)
results_ok_oligo_early = results_oligo_early[results_oligo_early$padjusted < 0.05,]
results_up_oligo_early = results_ok_oligo_early[results_ok_oligo_early$logFC > 0,]
results_down_oligo_early = results_ok_oligo_early[results_ok_oligo_early$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/logFC_barplot_results_oligo_early.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_oligo_early$logFC, decreasing = T))
dev.off()

res2_oligo = results_oligo_early %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene)

ranks_oligo = deframe(res2_oligo)
summary(ranks_oligo)

fgseaRes_oligo = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_oligo)
dim(fgseaRes_oligo)
colnames(fgseaRes_oligo)
table(fgseaRes_oligo$padj < 0.05)

fgseaResTidy_oligo = fgseaRes_oligo %>%
  as_tibble() %>%
  arrange(padj)

dim(fgseaResTidy_oligo)
fgseaResTidy_oligo


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/Hallmark_pathways_GSEA_oligo.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_oligo[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq oligo early---------------------------------------------------------------


genes_oligo=as.integer(p.adjust(results_oligo_early$pvalue[results_oligo_early$logFC!=0],method = "BH")<.05)

names(genes_oligo) = rownames(results_oligo_early[results_oligo_early$logFC!=0,])
table(genes_oligo)

pwf_oligo=nullp(genes_oligo,"hg19","geneSymbol",bias.data = NULL)

goResults_oligo = goseq(pwf_oligo, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_oligo)  #6099 funciones evaluadas (solo BP)
table(goResults_oligo$ontology)
table(goResults_oligo$over_represented_pvalue <0.05 |goResults_oligo$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/GO_results_oligo.jpeg", width=6, height=4, units="in", res=300)
goResults_oligo %>% 
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

#---------------------------------------------fGSEA oligoPC early---------------------------------------------------------------

table(results_oligoPC_early$padjusted < 0.05)
results_ok_oligoPC_early = results_oligoPC_early[results_oligoPC_early$padjusted < 0.05,]
results_up_oligoPC_early = results_ok_oligoPC_early[results_ok_oligoPC_early$logFC > 0,]
results_down_oligoPC_early = results_ok_oligoPC_early[results_ok_oligoPC_early$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/logFC_barplot_results_oligoPC_early.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_oligoPC_early$logFC, decreasing = T))
dev.off()

res2_oligoPC = results_oligoPC_early %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene)

ranks_oligoPC = deframe(res2_oligoPC)
summary(ranks_oligoPC)

fgseaRes_oligoPC = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_oligoPC)
dim(fgseaRes_oligoPC)
colnames(fgseaRes_oligoPC)
table(fgseaRes_oligoPC$padj < 0.05)

fgseaResTidy_oligoPC = fgseaRes_oligoPC %>%
  as_tibble() %>%
  arrange(padj)

dim(fgseaResTidy_oligoPC)
fgseaResTidy_oligoPC


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/Hallmark_pathways_GSEA_oligoPC.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_oligoPC[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq oligoPC early---------------------------------------------------------------


genes_oligoPC=as.integer(p.adjust(results_oligoPC_early$pvalue[results_oligoPC_early$logFC!=0],method = "BH")<.05)

names(genes_oligoPC) = rownames(results_oligoPC_early[results_oligoPC_early$logFC!=0,])
table(genes_oligoPC)

pwf_oligoPC=nullp(genes_oligoPC,"hg19","geneSymbol",bias.data = NULL)

goResults_oligoPC = goseq(pwf_oligoPC, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_oligoPC)  #6099 funciones evaluadas (solo BP)
table(goResults_oligoPC$ontology)
table(goResults_oligoPC$over_represented_pvalue <0.05 |goResults_oligoPC$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_early/GO_results_oligoPC.jpeg", width=6, height=4, units="in", res=300)
goResults_oligoPC %>% 
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
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------





