library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(Seurat)
library(msigdbr)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(presto)

setwd("D:/ZJU-FISH/BI/EC/Normal/")

# rm(list = ls())
sce <- readRDS("spe_nor_reset.rds")
# dge.celltype <- FindMarkers(sce, ident.1 = 'Memory CD4 T',ident.2 = 'Naive CD4 T', group.by = 'cell_type',logfc.threshold = 0,min.pct = 0)

DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce)
# sce <- subset(sce,ident=c('C0', 'C14'))
#dir.create("D:/ZJU-FISH/BI/EC/data/GSVA/kegg")
# setwd("D:/ZJU-FISH/BI/EC/data/GSVA/kegg")#选择基因集
# genesets <- msigdbr(species ="Homo sapiens", category = "C2")
# # unique(genesets$gs_subcat)
# genesets <- subset(genesets, gs_subcat == "CP:KEGG", select = c("gs_name","gene symbol")) %>% as.data.frame()
# genesets <- split(genesets$gene_symbol, genesets$gs_name)


genesets <- msigdbr(species = "Homo sapiens", category = "H")
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

Idents(sce) <- Idents(object = sce)
expr <- AverageExpression(sce, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
head(expr)

# gsva默认开启全部线程计算
gsva.res <- gsva(expr,genesets, method="gsva")
saveRDS(gsva.res,"gsva.spe.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names =F)
write.csv(gsva.df, "gsva_spe.csv", row.names = F)

pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row", cluster_cols = F)
