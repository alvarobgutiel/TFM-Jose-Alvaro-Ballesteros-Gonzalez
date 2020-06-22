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
library(EnsDb.Hsapiens.v79)


#-----------------------------------------------Loading data-----------------------------------------------

#Set working directory to the path of current file
dir_name = dirname(rstudioapi::getActiveDocumentContext()$path)
dir_real = strsplit(dir_name,split ="/scripts")[[1]]
setwd(dir_real)
dir.create("./figures/8.Gene Set Enrichment Analysis_20")

sce.hvg = readRDS("results/sce_diff_expr.rds")
results = readRDS("results/mostDE_list_20.rds")

#---------------------------------------------fGSEA---------------------------------------------------------------

table(results$padjusted < 0.05)
results_ok = results[results$padjusted < 0.05,]
results_up = results_ok[results_ok$logFC > 0,]
results_down = results_ok[results_ok$logFC < 0,]

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20/logFC_barplot.jpeg", width=6, height=4, units="in", res=300)
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

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20/Hallmark_pathways_GSEA_padj.jpeg", width=15, height=8, units="in", res=300)
ggplot(fgseaResTidy[1:25,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal()
dev.off()

#-----------------------------------------------goseq---------------------------------------------------------------

genes=as.integer(p.adjust(results$pvalue[results$logFC!=0],method = "BH")<.05)

names(genes) = rownames(results[results$logFC!=0,])
table(genes)

pwf=nullp(genes,"hg19","geneSymbol",bias.data = NULL)

goResults = goseq(pwf, "hg19","geneSymbol", test.cats=c("GO:BP"))
dim(goResults)  #6099 funciones evaluadas (solo BP)
table(goResults$ontology)
table(goResults$over_represented_pvalue <0.05 |goResults$under_represented_pvalue <0.05)

jpeg(file="./figures/8.Gene Set Enrichment Analysis_20/GO_results.jpeg", width=6, height=4, units="in", res=300)
goResults %>%
  top_n(15, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
  scale_fill_viridis_d(option="B")
dev.off()


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------





