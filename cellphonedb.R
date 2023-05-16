library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)

setwd("D:/ZJU-FISH/BI/EC/data/")
mysce = readRDS("des_escc.rds")
write.table(as.matrix(mysce@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(mysce@meta.data), mysce@meta.data[,'seurat_clusters', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA

write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)