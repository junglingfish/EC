library(SingleCellExperiment)
library(SCENIC)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(pheatmap)
library(foreach)
library(GENIE3)#推断基因共表达网络
library(SCopeLoomR)
library(loomR)
library(SeuratData)
library(SeuratDisk)

# setwd("D:/ZJU-FISH/BI/EC/data/")
# sceraw = readRDS('spe_nor_reset.rds')
# sce <- as.SingleCellExperiment(sce)
setwd("D:/ZJU-FISH/BI/EC/SCENIC_nor/")
# sce <- as.loom(sceraw, filename = "sceloom_nor.loom", verbose = FALSE)


# #下载
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather","https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# download.file(dbFiles[1], destfile = basename(dbFiles[1]), method = 'wget')
# download.file(dbFiles[1], destfile = basename(dbFiles[2]), method = 'wget')

#loom
# loomPath <- system.file(package="SCENIC", "sceloom.loom")
loom <- open_loom('D:/ZJU-FISH/BI/EC/SCENIC_nor/sceloom_nor.loom')
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=8) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# cd = sce[, Idents(sce) %in% c( "C0" , "C10")]
# exprMat <- counts(sce)
# cellInfo <- colData(sce)

#共表达网络
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
# exprMat_filtered_log <- log2(exprMat_filtered+1)
# runGenie3(exprMat_filtered, scenicOptions, nParts = 20)
runGenie3(exprMat_filtered, scenicOptions)
save(exprMat_filtered, file = "exprMat_filtered.rds")

#loom
loom <- open_loom('D:/ZJU-FISH/BI/EC/SCENIC_nor/sceloom_nor.loom')
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=1) 
saveRDS(scenicOptions, file="int/scenicOptions1.Rds") 

#
# setwd("/root/BI/EC/SCENIC/int/")
scenicOptions <- readRDS("int/scenicOptions1.Rds")
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

setwd("D:/ZJU-FISH/BI/EC/data/")
pbmc <- readRDS('spe_nor_reset.rds')
current.cluster.ids <- c(0:4)
new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4")
pbmc@meta.data$celltype <- plyr::mapvalues(x = pbmc@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
head(pbmc@meta.data)

##导入原始regulonAUC矩阵，也就是运行完runSCENIC_3_scoreCells后产生的AUC矩阵，可以查看每个GRN在每个细胞中的AUC活性打分
setwd("D:/ZJU-FISH/BI/EC/SCENIC/")
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- data.frame(t(AUCmatrix@assays@data@listData$AUC), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC) #把(替换成_
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC) #把)去掉，最后如把 KLF3_extended (79g) 替换成KLF3_extended_79g
colnames(AUCmatrix) <- RegulonName_AUC
pbmcauc <- AddMetaData(pbmc, AUCmatrix) #把AUC矩阵添加到pbmc的metadata信息中
pbmcauc@assays$integrated <- NULL
saveRDS(pbmcauc,'pbmcauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
pbmcbin <- AddMetaData(pbmc, BINmatrix)
pbmcbin@assays$integrated <- NULL
saveRDS(pbmcbin, 'pbmcbin.rds')

##利用Seurat可视化AUC
dir.create('scenic_seurat')
#FeaturePlot
GRNs <-intersect(colnames(AUCmatrix),colnames(BINmatrix))
for(i in 1:length(GRNs)){
  p1 = FeaturePlot(pbmcauc, features=GRNs[i], label=T, reduction = 'umap')
  p2 = FeaturePlot(pbmcbin, features=GRNs[i], label=T, reduction = 'umap')
  p3 = DimPlot(pbmc, reduction = 'umap', group.by = "celltype", label=T)
  plotc = p1|p2|p3
  ggsave(paste("scenic_seurat/",GRNs[i],".png",sep=""), plotc, width=14 ,height=4)
}

library(AUCell)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)

