# March 2020

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
library(EnsDb.Hsapiens.v79)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/8.Gene Set Enrichment Analysis_celltypesAD")

sce.hvg = readRDS("results/sce_diff_expr.rds")
results_astro = readRDS("results/mostDE_list_20_astro.rds")
results_endo = readRDS("results/mostDE_list_20_endo.rds")
results_micro = readRDS("results/mostDE_list_20_micro.rds")
results_oligo = readRDS("results/mostDE_list_20_oligo.rds")
results_oligoPC = readRDS("results/mostDE_list_20_oligoPC.rds")
results_neu = readRDS("results/mostDE_list_20_neu.rds")

pathways.hallmark = gmtPathways("data/msigdb.v7.1.symbols.gmt")

#---------------------------------------------fGSEA astro---------------------------------------------------------------
table(results_astro$padjusted < 0.05)
results_ok_astro = results_astro[results_astro$padjusted < 0.05,]
results_up_astro = results_ok_astro[results_ok_astro$logFC > 0,]
results_down_astro = results_ok_astro[results_ok_astro$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/logFC_barplot_astro.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_astro$logFC, decreasing = T))
dev.off()

res2_astro = results_astro %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_astro = deframe(res2_astro)
summary(ranks_astro)

fgseaRes_astro = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_astro, nperm=1000)
dim(fgseaRes_astro)
colnames(fgseaRes_astro)
table(fgseaRes_astro$padj < 0.05)

fgseaResTidy_astro = fgseaRes_astro %>%
  as_tibble() %>%
  arrange(desc(NES))

dim(fgseaResTidy_astro)
fgseaResTidy_astro

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/Hallmark_pathways_GSEA_astro.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_astro[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Astrocytes") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq astro---------------------------------------------------------------

genes_astro=as.integer(p.adjust(results_astro$pvalue[results_astro$logFC!=0],method = "BH")<.05)

names(genes_astro) = rownames(results_astro[results_astro$logFC!=0,])
table(genes_astro)

pwf_astro=nullp(genes_astro,"hg19","geneSymbol",bias.data = NULL)

goResults_astro = goseq(pwf_astro, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_astro)  #6099 funciones evaluadas (solo BP)
table(goResults_astro$ontology)
table(goResults_astro$over_represented_pvalue <0.05 |goResults_astro$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/GO_results_astro.jpeg", width=6, height=4, units="in", res=300)
goResults_astro %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")+
  ggtitle("GO results for Astrocytes")
dev.off()

#---------------------------------------------fGSEA endo---------------------------------------------------------------

table(results_endo$padjusted < 0.05)
results_ok_endo = results_endo[results_endo$padjusted < 0.05,]
results_up_endo = results_ok_endo[results_ok_endo$logFC > 0,]
results_down_endo = results_ok_endo[results_ok_endo$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/logFC_barplot_endo.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_endo$logFC, decreasing = T))
dev.off()

res2_endo = results_endo %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_endo = deframe(res2_endo)
summary(ranks_endo)

fgseaRes_endo = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_endo, nperm=1000)
dim(fgseaRes_endo)
colnames(fgseaRes_endo)
table(fgseaRes_endo$padj < 0.05)

fgseaResTidy_endo = fgseaRes_endo %>%
  as_tibble() %>%
  arrange(desc(NES))

dim(fgseaResTidy_endo)
fgseaResTidy_endo

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/Hallmark_pathways_GSEA_endo.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_endo[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Endothelial cells") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq endo---------------------------------------------------------------

genes_endo=as.integer(p.adjust(results_endo$pvalue[results_endo$logFC!=0],method = "BH")<.05)

names(genes_endo) = rownames(results_endo[results_endo$logFC!=0,])
table(genes_endo)

pwf_endo=nullp(genes_endo,"hg19","geneSymbol",bias.data = NULL)

goResults_endo = goseq(pwf_endo, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_endo)  #6099 funciones evaluadas (solo BP)
table(goResults_endo$ontology)
table(goResults_endo$over_represented_pvalue <0.05 |goResults_endo$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/GO_results_endo.jpeg", width=6, height=4, units="in", res=300)
goResults_endo %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")+
  ggtitle("GO results for Endothelial Cells")
dev.off()

#---------------------------------------------fGSEA micro---------------------------------------------------------------

table(results_micro$padjusted < 0.05)
results_ok_micro = results_micro[results_micro$padjusted < 0.05,]
results_up_micro = results_ok_micro[results_ok_micro$logFC > 0,]
results_down_micro = results_ok_micro[results_ok_micro$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/logFC_barplot_micro.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_micro$logFC, decreasing = T))
dev.off()

res2_micro = results_micro %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_micro = deframe(res2_micro)
summary(ranks_micro)

fgseaRes_micro = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_micro, nperm=1000)
dim(fgseaRes_micro)
colnames(fgseaRes_micro)
table(fgseaRes_micro$padj < 0.05)

fgseaResTidy_micro = fgseaRes_micro %>%
  as_tibble() %>%
  arrange(desc(NES))

dim(fgseaResTidy_micro)
fgseaResTidy_micro

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/Hallmark_pathways_GSEA_micro.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_micro[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Microglia") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq micro---------------------------------------------------------------

genes_micro=as.integer(p.adjust(results_micro$pvalue[results_micro$logFC!=0],method = "BH")<.05)

names(genes_micro) = rownames(results_micro[results_micro$logFC!=0,])
table(genes_micro)

pwf_micro=nullp(genes_micro,"hg19","geneSymbol",bias.data = NULL)

goResults_micro = goseq(pwf_micro, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_micro)  #6099 funciones evaluadas (solo BP)
table(goResults_micro$ontology)
table(goResults_micro$over_represented_pvalue <0.05 |goResults_micro$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/GO_results_micro.jpeg", width=6, height=4, units="in", res=300)
goResults_micro %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")+
  ggtitle("GO results for Microglia")
dev.off()

#---------------------------------------------fGSEA oligo---------------------------------------------------------------

table(results_oligo$padjusted < 0.05)
results_ok_oligo = results_oligo[results_oligo$padjusted < 0.05,]
results_up_oligo = results_ok_oligo[results_ok_oligo$logFC > 0,]
results_down_oligo = results_ok_oligo[results_ok_oligo$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/logFC_barplot_oligo.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_oligo$logFC, decreasing = T))
dev.off()

res2_oligo = results_oligo %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_oligo = deframe(res2_oligo)
summary(ranks_oligo)

fgseaRes_oligo = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_oligo, nperm=1000)
dim(fgseaRes_oligo)
colnames(fgseaRes_oligo)
table(fgseaRes_oligo$padj < 0.05)

fgseaResTidy_oligo = fgseaRes_oligo %>%
  as_tibble() %>%
  arrange(desc(NES))

dim(fgseaResTidy_oligo)
fgseaResTidy_oligo

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/Hallmark_pathways_GSEA_oligo1.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_oligo[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Oligodendrocytes") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq oligo---------------------------------------------------------------

genes_oligo=as.integer(p.adjust(results_oligo$pvalue[results_oligo$logFC!=0],method = "BH")<.05)

names(genes_oligo) = rownames(results_oligo[results_oligo$logFC!=0,])
table(genes_oligo)

pwf_oligo=nullp(genes_oligo,"hg19","geneSymbol",bias.data = NULL)

goResults_oligo = goseq(pwf_oligo, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_oligo)  #6099 funciones evaluadas (solo BP)
table(goResults_oligo$ontology)
table(goResults_oligo$over_represented_pvalue <0.05 |goResults_oligo$under_represented_pvalue <0.05)
#91 funciones significativas

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/GO_results_oligo.jpeg", width=6, height=4, units="in", res=300)
goResults_oligo %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")+
  ggtitle("GO results for Oligodendrocytes")
dev.off()

#---------------------------------------------fGSEA oligo PC---------------------------------------------------------------

table(results_oligoPC$padjusted < 0.05)
results_ok_oligoPC = results_oligoPC[results_oligoPC$padjusted < 0.05,]
results_up_oligoPC = results_ok_oligoPC[results_ok_oligoPC$logFC > 0,]
results_down_oligoPC = results_ok_oligoPC[results_ok_oligoPC$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/logFC_barplot_oligoPC.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_oligoPC$logFC, decreasing = T))
dev.off()

res2_oligoPC = results_oligoPC %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_oligoPC = deframe(res2_oligoPC)
summary(ranks_oligoPC)

fgseaRes_oligoPC = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_oligoPC, nperm=1000)
dim(fgseaRes_oligoPC)
colnames(fgseaRes_oligoPC)
table(fgseaRes_oligoPC$padj < 0.05)

fgseaResTidy_oligoPC = fgseaRes_oligoPC %>%
  as_tibble() %>%
  arrange(desc(NES))
dim(fgseaResTidy_oligoPC)
fgseaResTidy_oligoPC

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/Hallmark_pathways_GSEA_oligoPC.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_oligoPC[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Oligodendrocytes PC") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq oligo PC---------------------------------------------------------------

genes_oligoPC=as.integer(p.adjust(results_oligoPC$pvalue[results_oligoPC$logFC!=0],method = "BH")<.05)

names(genes_oligoPC) = rownames(results_oligoPC[results_oligoPC$logFC!=0,])
table(genes_oligoPC)

pwf_oligoPC=nullp(genes_oligoPC,"hg19","geneSymbol",bias.data = NULL)

goResults_oligoPC = goseq(pwf_oligoPC, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_oligoPC)  #6099 funciones evaluadas (solo BP)
table(goResults_oligoPC$ontology)
table(goResults_oligoPC$over_represented_pvalue <0.05 |goResults_oligoPC$under_represented_pvalue <0.05)
#91 funciones significativas

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/GO_results_oligoPC.jpeg", width=6, height=4, units="in", res=300)
goResults_oligoPC %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")+
  ggtitle("GO results for Oligodendrocytes PC")
dev.off()

#---------------------------------------------fGSEA neuron---------------------------------------------------------------

table(results_neu$padjusted < 0.05)
results_ok_neu = results_neu[results_neu$padjusted < 0.05,]
results_up_neu = results_ok_neu[results_ok_neu$logFC > 0,]
results_down_neu = results_ok_neu[results_ok_neu$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/logFC_barplot_neu.jpeg", width=6, height=4, units="in", res=300)
barplot(sort(results_neu$logFC, decreasing = T))
dev.off()

res2_neu = results_neu %>%
  dplyr::select(gene, logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene)

ranks_neu = deframe(res2_neu)
summary(ranks_neu)

fgseaRes_neu = fgsea(pathways=pathways.hallmark, minSize = 15, maxSize = 500, stats=ranks_neu, nperm=1000)
dim(fgseaRes_neu)
colnames(fgseaRes_neu)
table(fgseaRes_neu$padj < 0.05)

fgseaResTidy_neu = fgseaRes_neu %>%
  as_tibble() %>%
  arrange(desc(NES))

dim(fgseaResTidy_neu)
fgseaResTidy_neu

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/Hallmark_pathways_GSEA_neu.jpeg", width=15, height=10, units="in", res=300)
ggplot(fgseaResTidy_neu[1:40,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA Neurons") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq neuron---------------------------------------------------------------

genes_neu=as.integer(p.adjust(results_neu$pvalue[results_neu$logFC!=0],method = "BH")<.05)

names(genes_neu) = rownames(results_neu[results_neu$logFC!=0,])
table(genes_neu)

pwf_neu=nullp(genes_neu,"hg19","geneSymbol",bias.data = NULL)

goResults_neu = goseq(pwf_neu, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults_neu)  #6099 funciones evaluadas (solo BP)
table(goResults_neu$ontology)
table(goResults_neu$over_represented_pvalue <0.05 |goResults_neu$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_celltypesAD/GO_results_neu.jpeg", width=8, height=4, units="in", res=300)
goResults_neu %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")+
  ggtitle("GO results for Neurons")
dev.off()

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------





