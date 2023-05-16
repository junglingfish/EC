library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(tidyverse)
library(ROCR)
library(KernSmooth)
library(GEOquery)

# setwd("D:")
# file = read.csv('name.csv', sep = ',')
# name = subset(file, select = c('name'))
# name = as.vector(unlist(name[1]))
# # y = name[1, 1]
# # y
# 
# setwd("D:")
# file = read.csv('name.csv', sep = ',')
# p = subset(file, select = c('name'))
# p = as.vector(unlist(p))

setwd("D:/ZJU-FISH/BI/EC/data/GSE53624_RAW/")

# sce <- read.table(file="GSM1297076_US10313827_253831410058_S01_GE2_107_Sep09_1_1a.txt", header = TRUE, sep = '\t', quote = "", fill = T, comment.char = "!")
# head(sce)
# gse <- getGEO('gse53624', GSEMatrix = T, AnnotGPL = F)
# saveRDS(gse, file = "D:/ZJU-FISH/BI/EC/data/gse53624.rds")
# load('gse53624.Rdata')

sce <- read.table(file="GSM1297076_US10313827_253831410058_S01_GE2_107_Sep09_1_1a.txt")

mysce <- CreateSeuratObject(counts = sce, project = "escc", min.cells = 3, min.features = 200)
# head(mysce)
##mysce
mysce[["percent.mt"]] <- PercentageFeatureSet(mysce, pattern = "^MT-")
mysce <- subset(mysce, subset = nFeature_RNA < 2500 & nFeature_RNA > 400 & percent.mt < 10)


#双细胞去除
mysce <- NormalizeData(mysce)
mysce <- FindVariableFeatures(mysce, selection.method = "vst", nfeatures = 2000)
mysce <- ScaleData(mysce)
mysce <- RunPCA(mysce, verbose = FALSE)
mysce <- RunUMAP(mysce, dims = 1:15, verbose = FALSE)
mysce <- FindNeighbors(mysce, reduction = "pca", dims = 1:15)
mysce <- FindClusters(mysce, resolution = 0.5)
sweep.res.list_escc <- paramSweep_v3(mysce, PCs = 1:10, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
#print("FINISH!")
sweep.stats <- summarizeSweep(sweep.res.list_escc, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点

pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
#DoubletRate = ncol(pbmc)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
#DoubletRate = 0.05    # 直接查表，10000细胞对应的doublets rate是～7.6%
DoubletRate = ncol(mysce)*8*1e-6 #更通用
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞
homotypic.prop <- modelHomotypic(mysce$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi <- round(DoubletRate*ncol(mysce))
# 使用同源双细胞比例对计算的双细胞比例进行校正
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
mysce <- doubletFinder_v3(mysce, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

# #修改细胞名
# mysce <- RenameCells(mysce, new.names = name)

mysce <- NormalizeData(mysce, normalization.method = "LogNormalize", scale.factor = 1e4) 
mysce <- FindVariableFeatures(mysce, selection.method = 'vst', nfeatures = 2000)
mysce <- ScaleData(mysce, vars.to.regress = "percent.mt")
mysce <- RunPCA(mysce, features = VariableFeatures(object = mysce))

# mysce <- JackStraw(mysce, num.replicate = 100)
# mysce <- ScoreJackStraw(mysce, dims = 1:20)
# ElbowPlot(mysce)

mysce <- RunUMAP(mysce, reduction = "pca", dims = 1:15)
# mysce <- RunTSNE(mysce, reduction = "pca", dims = 1:15)
mysce <- FindNeighbors(mysce, dims = 1:15)
mysce <- FindClusters(mysce, resolution = 0.8)

p1 <- DimPlot(mysce, reduction = "umap")
p2 <- DimPlot(mysce, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

markers <- FindAllMarkers(mysce, only.pos = TRUE)
write.table(markers, file = 'Bcell_genes.csv', sep = ',')

# #TOTAL
features_escc <- c("CD19", "MS4A1", "IGHM", "IGHG3", "CCR7", "CD27", "CD1C", "CD24", "NR4A1", "CD99", "CD38", "BCL6", "CD14")
DotPlot(mysce, features = features_escc) + RotatedAxis()
# #eff
# features_escc <- c("CX3CR1", "KLRG1", "TBX21", "FCGR3A", "FGFBP2", "CXCR3", "CXCR5", "SLC4A10", "KLRB1", "RORC", "RORA")
# DotPlot(mysce, features = features_escc) + RotatedAxis()
# #eff
# features_escc <- c("HAVCR2", "CXCL13", "TOX")
# DotPlot(mysce, features = features_escc) + RotatedAxis()
# #NKT
# features_escc <- c("SLAMF1", "SLAMF6", "TGFBR1", "VA24", "JA18")
# DotPlot(mysce, features = features_escc) + RotatedAxis()
# #GDT
# features_escc <- c("IFNG", "IL17A", "IL17F", "IL22")
# DotPlot(mysce, features = features_escc) + RotatedAxis()
# #
# features_escc <- c("PDCD1", "CTLA4", "HAVCR2", "TIGIT", "LAG3", "CXCL13")
# DotPlot(mysce, features = features_escc) + RotatedAxis()
# 
# #top10
# mysce.markers <- FindAllMarkers(mysce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# top10 <- mysce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# 
# table(Idents(cd8))
# prop.table(table(Idents(mysce)))
# 
# FeaturePlot(mysce, features = c("CD3E", "CD4", "CD8A", "CD69", "CTLA4", "ITGAE", "GNLY", "NKG7", "PDCD1"), min.cutoff = "q9")

# FeaturePlot(mysce, features = c("P1T", "CD4", "CD8A", "CD69", "CTLA4", "ITGAE", "GNLY", "NKG7", "PDCD1"), min.cutoff = "q9")

# new.cluster.ids <- c("0", "1", "2", "3", "4", "0", "6", "7", "8", "9", "9", "11", "12", "13", "14", "15", "16", "2", "0")
# names(new.cluster.ids) <- levels(mysce)
# mysce <- RenameIdents(mysce, new.cluster.ids)
# DimPlot(mysce, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 5, raster=FALSE)

allcolour=c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#9932CC","#7FB4B2")

new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4", "C0", "C5", "C6", "C7", "C8", "C8", "C3", "C9", "C10", "C11", "C12", "C13", "C2", "C0")
names(new.cluster.ids) <- levels(mysce)
mysce <- RenameIdents(mysce, new.cluster.ids)
DimPlot(mysce, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6, raster=FALSE, cols = allcolour)

FeaturePlot(mysce, features = c("CD3D", "CD4", "CD8A", "ITGAE", "GNLY", "NKG7", "PDCD1", "CTLA4", "HAVCR2"), min.cutoff = "q9")

FeaturePlot(mysce, features = c("FOXP3", "IL2RA", "GZMK", "CCR7", "SELL", "TCF7", "IFNG", "LAYN", "ENTPD1"), min.cutoff = "q9")

FeaturePlot(mysce, features = c("CD38", "CD69", "MKI67", "STMN1", "MCM7", "GZMB", "IFNG", "CXCL13", "IL6ST"), min.cutoff = "q9")

#分期
DimPlot(object = mysce, group.by = "Sample", reduction.use = "umap", do.return = T, pt.size = 1, cells.use = grep("HC_1", TB.integrated@cell.names, value = T), cols.use = "green4") 
#

#complexheatmap
# cd48 = GetAssayData(object = mysce, slot = "counts")
library(ComplexHeatmap)
cdc = GetAssayData(mysce, slot = "counts")
cdc <- log2(cdc + 1)
gene_features <- c('CCL3', "GZMB", 'CD8A', 'VCAM1', 'LAG3', 'RHOB', 'RGS2', 'HAVCR2', 'PLPP1', 'TNF', 'TNFRSF4', 'FOXP3', 'IL2RA', 'LAIR2', 'TNFRSF18', 'BATF', 'CD177', 'IL1R2', 'CCR8', 'MAGEH1', 'GZMK', 'TUBA4A', 'CRTAM', 'CCL4L2', 'ITM2C', 'CXCR4', 'DUSP2', 'CMC1', 'TRAT1', 'ENC1', 'S100A10', 'S100A6', 'CD69', 'PERP', 'CD52', 'FKBP11', 'HSPA8', 'BTG2', 'HSPA6', 'ZFP36L2', 'IL7R', 'LMNA', 'ANXA1', 'GPR183', 'VIM', 'CDKN1A', 'MYADM', 'ZFP36', 'TPT1', 'RPS18', 'CCR7', 'PABPC1', 'PASK', 'SELL', 'KLF2', 'RPS3A', 'SLC2A2', 'CD55', 'FTH1', 'RPL32', 'CXCL13', 'NMB', 'NR3C1', 'TOX2', 'FKBP5', 'GK', 'PTPN13', 'GNG4', 'TNFSF8', 'IL6ST', 'TYROBP', 'FCER1G', 'XCL1', 'AREG', 'XCL2', 'GNLY', 'TRDC', 'KLRC1', 'KRT81', 'IFITM3', 'CCR6', 'TMEM173', 'PBXIP1', 'UCP2', 'SPOCK2', 'RORA', 'CCR4', 'CD28', 'FCMR', 'FGFBP2', 'SPON2', 'NKG7', 'PRF1', 'GZMH', 'S1PR5', 'FCGR3A', 'KLRD1', 'PLEK', 'KLRF1', 'STMN1', 'TYMS', 'TUBA1B', 'HMGN2', 'TUBB', 'DUT', 'RRM2', 'MKI67', 'HMGB2', 'UBE2C', 'KIAA0101', 'HSPA1B', 'BAG3', 'SERPINH1', 'HSPA1A', 'DNAJB1', 'ZFAND2A', 'HSPD1', 'CACYBP', 'HSPE1', 'IFIT2', 'IFIT3', 'ISG15', 'IFIT1', 'MX1', 'OASL', 'MX2', 'IFI6', 'RSAD2', 'ISG20', 'STAT1', 'HOPX', 'KLRC3', 'KIR3DL2', 'KLRC2', 'SYNGR1', 'LITAF', 'KIR2DL4', 'AOAH', 'PIK3R1', 'LDLRAD4')
cluster_info <- Idents(mysce)
cdc <- as.matrix(cdc[gene_features, names(cluster_info)])
Heatmap(cdc,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info)

#seurat heatmap
allcolor <- c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#9932CC","#7FB4B2")
gene_features <- c('CCL3', "GZMB", 'CD8A', 'VCAM1', 'LAG3', 'RHOB', 'RGS2', 'HAVCR2', 'PLPP1', 'TNF', 'TNFRSF4', 'FOXP3', 'IL2RA', 'LAIR2', 'TNFRSF18', 'BATF', 'CD177', 'IL1R2', 'CCR8', 'MAGEH1', 'GZMK', 'TUBA4A', 'CRTAM', 'CCL4L2', 'ITM2C', 'CXCR4', 'DUSP2', 'CMC1', 'TRAT1', 'ENC1', 'S100A10', 'S100A6', 'CD69', 'PERP', 'CD52', 'FKBP11', 'HSPA8', 'BTG2', 'HSPA6', 'ZFP36L2', 'IL7R', 'LMNA', 'ANXA1', 'GPR183', 'VIM', 'CDKN1A', 'MYADM', 'ZFP36', 'TPT1', 'RPS18', 'CCR7', 'PABPC1', 'PASK', 'SELL', 'KLF2', 'RPS3A', 'SLA2A2', 'CD55', 'FTH1', 'RPL32', 'CXCL13', 'NMB', 'NR3C1', 'TOX2', 'FKBP5', 'GK', 'PTPN13', 'GNG4', 'TNFSF8', 'IL6ST', 'TYROBP', 'FCER1G', 'XCL1', 'AREG', 'XCL2', 'GNLY', 'TRDC', 'KLRC1', 'KRT81', 'IFITM3', 'CCR6', 'TMEM173', 'PBXIP1', 'UCP2', 'SPOCK2', 'RORA', 'CCR4', 'CD28', 'FCMR', 'FGFBP2', 'SPON2', 'NKG7', 'PRF1', 'GZMH', 'S1PR5', 'FCGR3A', 'KLRD1', 'PLEK', 'KLRF1', 'STMN1', 'TYMS', 'TUBA1B', 'HMGN2', 'TUBB', 'DUT', 'RRM2', 'MKI67', 'HMGB2', 'UBE2C', 'KIAA0101', 'HSPA1B', 'BAG3', 'SERPINH1', 'HSPA1A', 'DNAJB1', 'ZFAND2A', 'HSPD1', 'CACYBP', 'HSPE1', 'IFIT2', 'IFIT3', 'ISG15', 'IFIT1', 'MX1', 'OASL', 'MX2', 'IFI6', 'RSAD2', 'ISG20', 'STAT1', 'HOPX', 'KLRC3', 'KIR3DL2', 'KLRC2', 'SYNGR1', 'LITAF', 'KIR2DL4', 'AOAH', 'PIK3R1', 'LDLRAD4')
DoHeatmap(mysce, features = gene_features, group.colors = allcolor)


features_c <- c("CD2", "CD3D", "CD3E", "CD3G", "ITGAE", "PDCD1", "CTLA4", "CD19", "CD79A", "MS4A1", "JCHAIN", "MZB1", "TRDC", "KLRD1", "GNLY", "NKG7")
p <- DoHeatmap(mysce, features = features_c, label = F, group.colors = 'white')
p

features_escc <- c("CD69", "ITGAE", "CD4", "CD8A", "PDCD1", "CTLA4", "CCR7", "SELL", "TCF7", "LEF1", "FOXP3", "GZMK", "CD3D", "GNLY", "NKG7", "KLRD1", "ENTPD1", "CD38")
DotPlot(mysce, features = features_escc) + RotatedAxis()
