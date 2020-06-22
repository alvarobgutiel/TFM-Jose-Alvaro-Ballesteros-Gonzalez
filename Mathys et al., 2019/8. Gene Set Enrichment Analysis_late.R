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
dir.create("./figures/8.Gene Set Enrichment Analysis_20_cell_type_late")

sce.hvg = readRDS("results/sce_diff_expr.rds")
results_astro_late = readRDS("results/mostDE_list_20_astro_late.rds")
results_micro_late = readRDS("results/mostDE_list_20_micro_late.rds")
results_endo_late = readRDS("results/mostDE_list_20_endo_late.rds")
results_neu_late = readRDS("results/mostDE_list_20_neu_late.rds")
results_oligo_late = readRDS("results/mostDE_list_20_oligo_late.rds")
results_oligoPC_late = readRDS("results/mostDE_list_20_oligoPC_late.rds")

#---------------------------------------------fGSEA astro late---------------------------------------------------------------

table(results_astro_late$padjusted < 0.05)
results_ok_astro_late = results_astro_late[results_astro_late$padjusted < 0.05,]
results_up_astro_late = results_ok_astro_late[results_ok_astro_late$logFC > 0,]
results_down_astro_late = results_ok_astro_late[results_ok_astro_late$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/logFC_barplot_results_astro_late.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_astro_late$logFC, decreasing = T))
dev.off()

res2_astro = results_astro_late %>% 
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


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/Hallmark_pathways_GSEA_astro.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_astro[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq astro late---------------------------------------------------------------


genes_astro=as.integer(p.adjust(results_astro_late$pvalue[results_astro_late$logFC!=0],method = "BH")<.05)

names(genes_astro) = rownames(results_astro_late[results_astro_late$logFC!=0,])
table(genes_astro)

pwf_astro=nullp(genes_astro,"hg19","geneSymbol",bias.data = NULL)

goResults_astro = goseq(pwf_astro, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_astro)  #6099 funciones evaluadas (solo BP)
table(goResults_astro$ontology)
table(goResults_astro$over_represented_pvalue <0.05 |goResults_astro$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/GO_results_astro.jpeg", width=6, height=4, units="in", res=300)
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

#---------------------------------------------fGSEA micro late---------------------------------------------------------------

table(results_micro_late$padjusted < 0.05)
results_ok_micro_late = results_micro_late[results_micro_late$padjusted < 0.05,]
results_up_micro_late = results_ok_micro_late[results_ok_micro_late$logFC > 0,]
results_down_micro_late = results_ok_micro_late[results_ok_micro_late$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/logFC_barplot_results_micro_late.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_micro_late$logFC, decreasing = T))
dev.off()

res2_micro = results_micro_late %>% 
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


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/Hallmark_pathways_GSEA_micro.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_micro[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +  scale_fill_manual(values=c("#00BFC4"))+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq micro late---------------------------------------------------------------


genes_micro=as.integer(p.adjust(results_micro_late$pvalue[results_micro_late$logFC!=0],method = "BH")<.05)

names(genes_micro) = rownames(results_micro_late[results_micro_late$logFC!=0,])
table(genes_micro)

pwf_micro=nullp(genes_micro,"hg19","geneSymbol",bias.data = NULL)

goResults_micro = goseq(pwf_micro, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_micro)  #6099 funciones evaluadas (solo BP)
table(goResults_micro$ontology)
table(goResults_micro$over_represented_pvalue <0.05 |goResults_micro$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/GO_results_micro.jpeg", width=6, height=4, units="in", res=300)
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

#---------------------------------------------fGSEA endo late---------------------------------------------------------------

table(results_endo_late$padjusted < 0.05)
results_ok_endo_late = results_endo_late[results_endo_late$padjusted < 0.05,]
results_up_endo_late = results_ok_endo_late[results_ok_endo_late$logFC > 0,]
results_down_endo_late = results_ok_endo_late[results_ok_endo_late$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/logFC_barplot_results_endo_late.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_endo_late$logFC, decreasing = T))
dev.off()

res2_endo = results_endo_late %>% 
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


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/Hallmark_pathways_GSEA_endo.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_endo[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq endo late---------------------------------------------------------------


genes_endo=as.integer(p.adjust(results_endo_late$pvalue[results_endo_late$logFC!=0],method = "BH")<.05)

names(genes_endo) = rownames(results_endo_late[results_endo_late$logFC!=0,])
table(genes_endo)

pwf_endo=nullp(genes_endo,"hg19","geneSymbol",bias.data = NULL)

goResults_endo = goseq(pwf_endo, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_endo)  #6099 funciones evaluadas (solo BP)
table(goResults_endo$ontology)
table(goResults_endo$over_represented_pvalue <0.05 |goResults_endo$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/GO_results_endo.jpeg", width=6, height=4, units="in", res=300)
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

#---------------------------------------------fGSEA neu late---------------------------------------------------------------

table(results_neu_late$padjusted < 0.05)
results_ok_neu_late = results_neu_late[results_neu_late$padjusted < 0.05,]
results_up_neu_late = results_ok_neu_late[results_ok_neu_late$logFC > 0,]
results_down_neu_late = results_ok_neu_late[results_ok_neu_late$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/logFC_barplot_results_neu_late.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_neu_late$logFC, decreasing = T))
dev.off()

res2_neu = results_neu_late %>% 
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
  arrange(desc(NES))

dim(fgseaResTidy_neu)
fgseaResTidy_neu


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/Hallmark_pathways_GSEA_neu.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_neu[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05))+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq neu late---------------------------------------------------------------


genes_neu=as.integer(p.adjust(results_neu_late$pvalue[results_neu_late$logFC!=0],method = "BH")<.05)

names(genes_neu) = rownames(results_neu_late[results_neu_late$logFC!=0,])
table(genes_neu)

pwf_neu=nullp(genes_neu,"hg19","geneSymbol",bias.data = NULL)

goResults_neu = goseq(pwf_neu, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_neu)  #6099 funciones evaluadas (solo BP)
table(goResults_neu$ontology)
table(goResults_neu$over_represented_pvalue <0.05 |goResults_neu$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/GO_results_neu.jpeg", width=6, height=4, units="in", res=300)
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

#---------------------------------------------fGSEA oligo late---------------------------------------------------------------

table(results_oligo_late$padjusted < 0.05)
results_ok_oligo_late = results_oligo_late[results_oligo_late$padjusted < 0.05,]
results_up_oligo_late = results_ok_oligo_late[results_ok_oligo_late$logFC > 0,]
results_down_oligo_late = results_ok_oligo_late[results_ok_oligo_late$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/logFC_barplot_results_oligo_late.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_oligo_late$logFC, decreasing = T))
dev.off()

res2_oligo = results_oligo_late %>% 
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


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/Hallmark_pathways_GSEA_oligo.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_oligo[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq oligo late---------------------------------------------------------------


genes_oligo=as.integer(p.adjust(results_oligo_late$pvalue[results_oligo_late$logFC!=0],method = "BH")<.05)

names(genes_oligo) = rownames(results_oligo_late[results_oligo_late$logFC!=0,])
table(genes_oligo)

pwf_oligo=nullp(genes_oligo,"hg19","geneSymbol",bias.data = NULL)

goResults_oligo = goseq(pwf_oligo, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_oligo)  #6099 funciones evaluadas (solo BP)
table(goResults_oligo$ontology)
table(goResults_oligo$over_represented_pvalue <0.05 |goResults_oligo$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/GO_results_oligo.jpeg", width=6, height=4, units="in", res=300)
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

#---------------------------------------------fGSEA oligoPC late---------------------------------------------------------------

table(results_oligoPC_late$padjusted < 0.05)
results_ok_oligoPC_late = results_oligoPC_late[results_oligoPC_late$padjusted < 0.05,]
results_up_oligoPC_late = results_ok_oligoPC_late[results_ok_oligoPC_late$logFC > 0,]
results_down_oligoPC_late = results_ok_oligoPC_late[results_ok_oligoPC_late$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/logFC_barplot_results_oligoPC_late.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_oligoPC_late$logFC, decreasing = T))
dev.off()

res2_oligoPC = results_oligoPC_late %>% 
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


jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/Hallmark_pathways_GSEA_oligoPC.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_oligoPC[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astro") + 
  theme_minimal()
dev.off()

#-----------------------------------------------goseq oligoPC late---------------------------------------------------------------


genes_oligoPC=as.integer(p.adjust(results_oligoPC_late$pvalue[results_oligoPC_late$logFC!=0],method = "BH")<.05)

names(genes_oligoPC) = rownames(results_oligoPC_late[results_oligoPC_late$logFC!=0,])
table(genes_oligoPC)

pwf_oligoPC=nullp(genes_oligoPC,"hg19","geneSymbol",bias.data = NULL)

goResults_oligoPC = goseq(pwf_oligoPC, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_oligoPC)  #6099 funciones evaluadas (solo BP)
table(goResults_oligoPC$ontology)
table(goResults_oligoPC$over_represented_pvalue <0.05 |goResults_oligoPC$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20_cell_type_late/GO_results_oligoPC.jpeg", width=6, height=4, units="in", res=300)
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





