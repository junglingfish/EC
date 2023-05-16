library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(tidyverse)
library(ROCR)
library(KernSmooth)

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

setwd("D:/ZJU-FISH/BI/EC/data/")

sce <- read.table(file="UMI_matrix_Bcell.txt")
#head(sce)
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

sce_cd4 = mysce[, Idents(mysce) %in% c( "C1" , "C4", "C5", "C6", "C8", "C11" )]
saveRDS(sce_cd4, file = "cd4.rds")

# y = Idents(sce_cd4)
# write.table(y, file = "cd4.csv", sep = ',')

sce_cd8 = mysce[, Idents(mysce) %in% c( "C0" , "C2", "C10", "C12" )]
saveRDS(sce_cd8, file = "cd8.rds")

# y = Idents(sce_cd8)
# write.table(y, file = "cd8.csv", sep = ',')

#病人基因表达量
cellname = colnames(mysce)
scemat = GetAssayData(mysce, slot = "counts")
cell1 = cellname[1:347]
scemat1 = scemat[][1:347]


library(limma)
library(edgeR)
library(statmod)
#差异基因
cd010 = mysce[, Idents(mysce) %in% c( "C0" , "C1")]
cell = GetAssayData(object = cd010, slot = "counts")
dgelist <- DGEList(counts = cell, group = Idents(cd010))
keep <- rowSums(cpm(dgelist) > 1 ) >= 2 #过滤
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')  #TMM标准化
group <- rep(c('C0', 'C1'), each = 12850, length.out = 19148) 
design <- model.matrix(~group)    #构建分组矩阵
dge <- estimateDisp(dgelist_norm, design, robust = TRUE) #估算离散值

#差异分析
fit <- glmFit(dge, design, robust = TRUE) #拟合模型
lrt <- glmLRT(fit)   #统计检验
topTags(lrt)

write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'deg.csv', quote = FALSE) #输出主要结果

#volcano map
diff_stat <- read.csv("c0_c1_deg.csv", header = TRUE, row.names = 1)
diff_stat[which(diff_stat$FDR < 0.05 & diff_stat$logFC >= 1),'diff'] <- 'up' #上调趋势筛选
diff_stat[which(diff_stat$FDR < 0.05 & diff_stat$logFC <= -1),'diff'] <- 'down' #下调趋势筛选
diff_stat[!(diff_stat$diff %in% c('up', 'down')),'diff'] <- 'no'
p1 <- ggplot(diff_stat, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = diff), size = 1) +
  scale_colour_manual(limits = c('up', 'down', 'no'), values = c('blue', 'red', 'gray40'), labels = c('Enriched OTUs', 'Depleted OTUs', 'No diff OTUs')) +
  labs(x = 'log2 Fold Change', y = '-log10 FDR p-value')
p1

#cd0 cd10分群
cd = readRDS('des_escc.rds')
cdall = cd[, Idents(cd) %in% c( "C0" , "C10")]

cdall <- NormalizeData(cdall, normalization.method = "LogNormalize", scale.factor = 1e4) 
cdall <- FindVariableFeatures(cdall, selection.method = 'vst', nfeatures = 2000)
cdall <- ScaleData(cdall, vars.to.regress = "percent.mt")
cdall <- RunPCA(cdall, features = VariableFeatures(object = cdall))

# cdall <- JackStraw(cdall, num.replicate = 100)
# cdall <- ScoreJackStraw(cdall, dims = 1:20)
# ElbowPlot(cdall)

#15
cdall <- RunUMAP(cdall, reduction = "pca", dims = 1:15)
# cdall <- RunTSNE(cdall, reduction = "pca", dims = 1:15)
cdall <- FindNeighbors(cdall, dims = 1:15)
cdall <- FindClusters(cdall, resolution = 0.2)

p1 <- DimPlot(cdall, reduction = "umap")
p2 <- DimPlot(cdall, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4")
names(new.cluster.ids) <- levels(cdall)
cdall <- RenameIdents(cdall, new.cluster.ids)
DimPlot(cdall, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6, raster=FALSE)

markers <- FindAllMarkers(cdall, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file = 'c0_c10_all.csv', sep = ',')

###CD8分群
cd = readRDS('des_escc.rds')
cd8 = cd[, Idents(cd) %in% c( "C0" , "C2", "C3", "C9", "C10", "C12")]

cd8 <- NormalizeData(cd8, normalization.method = "LogNormalize", scale.factor = 1e4) 
cd8 <- FindVariableFeatures(cd8, selection.method = 'vst', nfeatures = 2000)
cd8 <- ScaleData(cd8, vars.to.regress = "percent.mt")
cd8 <- RunPCA(cd8, features = VariableFeatures(object = cd8))

# cd8 <- JackStraw(cd8, num.replicate = 100)
# cd8 <- ScoreJackStraw(cd8, dims = 1:20)
# ElbowPlot(cd8)

#11
cd8 <- RunUMAP(cd8, reduction = "pca", dims = 1:11)
# cd8 <- RunTSNE(cd8, reduction = "pca", dims = 1:15)
cd8 <- FindNeighbors(cd8, dims = 1:11)
cd8 <- FindClusters(cd8, resolution = 0.15)

# #12
# cd8 <- RunUMAP(cd8, reduction = "pca", dims = 1:12)
# # cd8 <- RunTSNE(cd8, reduction = "pca", dims = 1:15)
# cd8 <- FindNeighbors(cd8, dims = 1:12)
# cd8 <- FindClusters(cd8, resolution = 0.3)

p1 <- DimPlot(cd8, reduction = "umap")
p2 <- DimPlot(cd8, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#
features_escc <- c('CD69', 'ITGAE', 'ITGA1', 'PDCD1', 'CTLA4', 'CD38', 'ENTPD1', 'MKI67')
DotPlot(cd8, features = features_escc) + RotatedAxis()

#
features_escc <- c("CD69", "ITGAE", "S1PR1", "KLF2", "CCR7", "ZNF683", "SELL", "EOMES", "PRDM1")
DotPlot(cd8, features = features_escc) + RotatedAxis()

FeaturePlot(cd8, features = c("CD69", "ITGAE", "S1PR1", "KLF2", "CCR7", "ZNF683", "SELL", "EOMES", "PRDM1"), min.cutoff = "q9")

#
features_escc <- c("GZMK","PDCD1","HAVCR2","CD101","CD8A","CTLA4","GZMB","ITGA1","CXCR6")
DotPlot(cd8, features = features_escc) + RotatedAxis()

FeaturePlot(cd8, features = c("GZMK","PDCD1","HAVCR2","CD101","CD8A","CTLA4","GZMB","ITGA1","CXCR6"), min.cutoff = "q9")

new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6")
names(new.cluster.ids) <- levels(cd8)
cd8 <- RenameIdents(cd8, new.cluster.ids)
DimPlot(cd8, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6, raster=FALSE)

#cd8 label clusters
tmem = c(0)
tex = c(1)
trm = c(2)
trmpro = c(3)
temra = c(4)
tnaive = c(5)
tcyto = c(6)

current.cluster.ids <- c(tmem, tex, trm, trmpro, temra, tnaive, tcyto)

new.cluster.ids <- c(rep("Tmem",length(tmem)), rep("Tex",length(tex)), rep("TRM",length(trm)), rep("TRM_pro",length(trmpro)), rep("Temra",length(temra)), rep("Tnaive",length(tnaive)), rep("Tcyto",length(tcyto)))

mysce@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(mysce@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

saveRDS(mysce, file = 'cd8_labeled.rds')

#cd8 heatmap
allcolor <- c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080")
gene_features <- c('IL7R', "LMNA", 'MT1X', 'GPR183', 'MT2A', 'ANXA1', 'MYADM', 'FOS', 'TUBA4A', 'FOSB', 'CXCL13', 'KRT86', 'SOX4', 'TIGIT', 'PLPP1', 'ACP5', 'LAYN', 'CXCR6', 'ENTPD1', 'CTLA4', 'HSPA6', 'SERPINH1', 'RHOB', 'HSPA1B', 'BAG3', 'HSPB1', 'CCL3', 'ID3', 'HMOX1', 'DNAJB1', 'STMN1', 'TYMS', 'TUBB', 'TUBA1B', 'HMGN2', 'RRM2', 'DUT', 'UBE2C', 'MKI67', 'KIAA0101', 'FGFBP2', 'TYROBP', 'SPON2', 'FCGR3A', 'S1PR5', 'PLAC8', 'KLRF1', 'PLEK', 'KLRD1', 'NKG7', 'CCR7', 'SELL', 'LTB', 'ISG15', 'MX1', 'MX2', 'SAT1', 'IFI6', 'TNFRSF4', 'IFIT3', 'CRTAM', 'CHI3L2', 'CAV1', 'GZMK', 'XCL2', 'NMB', 'CD27', 'ZNF331', 'HES4', 'DUSP4')
DoHeatmap(cd8, features = gene_features, group.colors = allcolor)+scale_fill_gradientn(colors = c("#1A5B86","white","#921422"))#设置热图颜色  

markers <- FindAllMarkers(cd8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file = 'cd8_all.csv', sep = ',')

#result
features <- c('CD69', 'ITGAE', 'ITGA1', 'CCR7', 'SELL', 'S1PR1')
# DotPlot(mysce, features = features) + RotatedAxis()
DotPlot(cd8, features = features) + RotatedAxis()

saveRDS(cd8, file = "cd8.rds")

# #violin
# VlnPlot(cd8, features = c('HAVCR2', 'PDCD1', 'LAG3', 'TIGIT', 'CTLA4', 'LAYN', 'ENTPD1'))
library(MySeuratWrappers)  
#需要展示的基因  
markers <- c('HAVCR2', 'PDCD1', 'LAG3', 'TIGIT', 'CTLA4', 'LAYN', 'ENTPD1')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(cd8, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度

markers <- c('PRDM1', 'ZNF683', 'RUNX3', 'TBX21', 'EOMES', 'KLF2', 'TOX')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(cd8, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度

markers <- c('CD69', 'ITGAE', 'ITGA1', 'PDCD1', 'CTLA4', 'CD38', 'ENTPD1', 'MKI67')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(cd8, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度

markers <- c('LAG3', 'TIGIT', 'PDCD1', 'HAVCR2', 'CTLA4', 'LAYN', 'NKG2A')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(cd8, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度

markers <- c('PRDM1', 'ZNF683', 'RNUX3', 'TBX21', 'EOMES', 'KLF2', 'TOX', 'ID2', 'ID3', 'BATF', 'MAF', 'RBPJ', 'BHLHE40')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(cd8, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度


#point cd8
cd8 = readRDS('cd8.rds')
cd8s = cd8[, Idents(cd8) %in% c( "C1" , "C2", "C3", "C6")]
cd8s <- NormalizeData(cd8s, normalization.method = "LogNormalize", scale.factor = 1e4)
cd8s <- FindVariableFeatures(cd8s, selection.method = 'vst', nfeatures = 2000)
cd8s <- ScaleData(cd8s, vars.to.regress = "percent.mt")
cd8s <- RunPCA(cd8s, features = VariableFeatures(object = cd8s))

#11
cd8s <- RunUMAP(cd8s, reduction = "pca", dims = 1:10)
# cd8s <- RunTSNE(cd8s, reduction = "pca", dims = 1:15)
cd8s <- FindNeighbors(cd8s, dims = 1:10)
cd8s <- FindClusters(cd8s, resolution = 0.23)

p1 <- DimPlot(cd8s, reduction = "umap")
p2 <- DimPlot(cd8s, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#
features_escc <- c("CD69", "ITGAE", "CD38", "ENTPD1", "ITGA1", "PDCD1", "CD27", "CD28", "CTLA4", "CD101", 'MKI67', 'EOMES')
DotPlot(cd8, features = features_escc) + RotatedAxis()
DotPlot(cd8s, features = features_escc)+coord_flip()+theme_bw()+theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066','#336699','#66CC66','#FFCC33'))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

###根据cd8 4群
cd8 = readRDS('cd8.rds')
cd84 = cd8[, Idents(cd8) %in% c("C4")]

cd84 <- NormalizeData(cd84, normalization.method = "LogNormalize", scale.factor = 1e4)
cd84 <- FindVariableFeatures(cd84, selection.method = 'vst', nfeatures = 2000)
cd84 <- ScaleData(cd84, vars.to.regress = "percent.mt")
cd84 <- RunPCA(cd84, features = VariableFeatures(object = cd84))

#11
cd84 <- RunUMAP(cd84, reduction = "pca", dims = 1:11)
# cd84 <- RunTSNE(cd84, reduction = "pca", dims = 1:15)
cd84 <- FindNeighbors(cd84, dims = 1:11)
cd84 <- FindClusters(cd84, resolution = 0.15)

p1 <- DimPlot(cd84, reduction = "umap")
p2 <- DimPlot(cd84, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#
features_escc <- c("KLRG1")
DotPlot(cd84, features = features_escc) + RotatedAxis()

new.cluster.ids <- c("C0", "C1", "C1")
names(new.cluster.ids) <- levels(cd84)
cd84 <- RenameIdents(cd84, new.cluster.ids)
DimPlot(cd84, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, raster=FALSE)

# exp = FindConservedMarkers(cd84, ident.1 = 0,  ident.2 = 1, grouping.var = 'seurat_clusters', only.pos = TRUE,
#                            logfc.threshold = 0.25)
# # cd84_1 <- WhichCells(object = cd84, expression = KLRG1 > 0.5)
# cd84_1 <- subset(cd84, subset = KLRG1 >= 0.5)
# cd84_2 <- subset(cd84, subset = KLRG1 < 0.5)

markers <- FindAllMarkers(cd84, only.pos = TRUE)
write.table(markers, file = 'cd84.csv', sep = ',')

# markers <- FindAllMarkers(cd84_2, only.pos = TRUE)
# write.table(markers, file = 'cd84_2.csv', sep = ',')

###cd4
cd = readRDS('des_escc.rds')
cd4 = cd[, Idents(cd) %in% c( "C1" , "C4", "C5", "C6", "C8", "C11")]

cd4 <- NormalizeData(cd4, normalization.method = "LogNormalize", scale.factor = 1e4) 
cd4 <- FindVariableFeatures(cd4, selection.method = 'vst', nfeatures = 2000)
cd4 <- ScaleData(cd4, vars.to.regress = "percent.mt")
cd4 <- RunPCA(cd4, features = VariableFeatures(object = cd4))

# cd4 <- JackStraw(cd4, num.replicate = 100)
# cd4 <- ScoreJackStraw(cd4, dims = 1:20)
# ElbowPlot(cd4)

#11
cd4 <- RunUMAP(cd4, reduction = "pca", dims = 1:15)
# cd4 <- RunTSNE(cd4, reduction = "pca", dims = 1:15)
cd4 <- FindNeighbors(cd4, dims = 1:15)
cd4 <- FindClusters(cd4, resolution = 0.35)

# #12
# cd4 <- RunUMAP(cd4, reduction = "pca", dims = 1:12)
# # cd4 <- RunTSNE(cd4, reduction = "pca", dims = 1:15)
# cd4 <- FindNeighbors(cd4, dims = 1:12)
# cd4 <- FindClusters(cd4, resolution = 0.3)

p1 <- DimPlot(cd4, reduction = "umap")
p2 <- DimPlot(cd4, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

markers <- FindAllMarkers(cd4, only.pos = TRUE)
write.table(markers, file = 'cd4.csv', sep = ',')

#
features_escc <- c("SELL", "ITGAE", "S1PR1", "KLF2", "CCR7", "LEF1", "TCF7", "KLRG1", "CX3CR1", 'CTLA4', "IL23R", "IL17A")
DotPlot(cd4, features = features_escc) + RotatedAxis()

#
features_escc <- c("IL2RA", "FOXP3", "CCR8", "PDCD1", "HAVCR2", "IFNG", "GZMK", "CXCR5", "CXCR6", "GNLY", "IL10", "CXCL13")
DotPlot(cd4, features = features_escc) + RotatedAxis()

#CD4 heatmap
allcolor <- c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080")
gene_features <- c('TNFRSF4', 'FOXP3', 'IL2RA', 'TNFRSF9', 'TNFRSF18', 'CARD16', 'IFI6', 'IL32', 'LAYN', 'ANXA1', 'CCR7', 'FOS', 'KLF2', 'CD55', 'ZFP36L2', 'RPS12', 'S1PR1', 'MYADM', 'IL7R', 'GZMA', 'KLRB1', 'IL17A', 'CAPG', 'CTSH', 'MGAT4A', 'ZEB2', 'FURIN', 'PKIG', 'IL26', 'CXCL13', 'CCL4', 'G0S2', 'NMB', 'ALOX5AP', 'NR3C1', 'BHLHE40', 'IFNG', 'CD200', 'KLRB1', 'HSPA6', 'HSPA1B', 'SERPINH1', 'BAG3', 'HSPA1A', 'ZFAND2A', 'HSPB1', 'DNAJB1', 'DNAJB4', 'HSPD1', 'HLA-DRB1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB5', 'HLA-DPB1', 'CD74', 'SELL', 'SAMHD1', 'LEF1', 'KLF3', 'TOX2', 'GNG4', 'CHI3L2', 'CXCR5', 'TCF7', 'SMCO4', 'TOX', 'PTPN13', 'ICA1', 'POU2AF1', 'CCL5', 'NKG7', 'CTSW', 'GZMK', 'YBX3', 'KLRD1', 'ZNF683', 'KLRC1', 'CRTAM', 'GNLY')
DoHeatmap(cd4, features = gene_features, group.colors = allcolor)+scale_fill_gradientn(colors = c("#1A5B86","white","#921422"))#设置热图颜色  


#cd4 label clusters
treg1 = c(0)
naive = c(1)
th17 = c(2)
th1 = c(3)
treg2 = c(4)
treg3 = c(5)
tfh = c(6)
tcm = c(7)

current.cluster.ids <- c(treg1, naive, th17, th1, treg2, treg3, tfh, tcm)

new.cluster.ids <- c(rep("Treg1",length(treg1)), rep("Tnaive",length(naive)), rep("Th17",length(th17)), rep("Th1_like",length(th1)), rep("Treg2",length(treg2)), rep("Treg3",length(treg3)), rep("Tfh",length(tfh)), rep("Tcm",length(tcm)))

mysce@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(mysce@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

saveRDS(mysce, file = 'cd4_labeled.rds')

### 提取特定基因
mysce <- readRDS('des_escc.rds')
select_cells <- WhichCells(mysce, slot = 'counts', expression = CD3D > 0 & CD8A > 0 & CD69 > 0 & ITGAE > 0)
select_obj <- subset(mysce, cells = select_cells)

select_obj <- NormalizeData(select_obj, normalization.method = "LogNormalize", scale.factor = 1e4) 
select_obj <- FindVariableFeatures(select_obj, selection.method = 'vst', nfeatures = 2000)
select_obj <- ScaleData(select_obj, vars.to.regress = "percent.mt")
select_obj <- RunPCA(select_obj, features = VariableFeatures(object = select_obj))

# select_obj <- JackStraw(select_obj, num.replicate = 100)
# select_obj <- ScoreJackStraw(select_obj, dims = 1:20)
# ElbowPlot(select_obj)

#11
select_obj <- RunUMAP(select_obj, reduction = "pca", dims = 1:15)
# select_obj <- RunTSNE(select_obj, reduction = "pca", dims = 1:15)
select_obj <- FindNeighbors(select_obj, dims = 1:15)
select_obj <- FindClusters(select_obj, resolution = 0.25)

# #12
# select_obj <- RunUMAP(select_obj, reduction = "pca", dims = 1:12)
# # select_obj <- RunTSNE(select_obj, reduction = "pca", dims = 1:15)
# select_obj <- FindNeighbors(select_obj, dims = 1:12)
# select_obj <- FindClusters(select_obj, resolution = 0.3)

p1 <- DimPlot(select_obj, reduction = "umap")
p2 <- DimPlot(select_obj, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#
features_escc <- c('CD69', 'ITGAE', 'ITGA1', 'PDCD1', 'CTLA4', 'EOMES', 'MKI67', 'CD38', 'ENTPD1', 'GZMB', 'GNLY', 'IFNG', 'CCL5', 'TNF', 'GZMA', 'CD101')
DotPlot(select_obj, features = features_escc) + RotatedAxis()

#
features_escc <- c('CD3E', 'CD8A', 'CXCR6')
DotPlot(select_obj, features = features_escc) + RotatedAxis()

#result
features <- c('CD69', 'ITGAE', 'ITGA1', 'CCR7', 'SELL', 'S1PR1')
# DotPlot(mysce, features = features) + RotatedAxis()
DotPlot(select_obj, features = features) + RotatedAxis()

# #violin
# VlnPlot(cd8, features = c('HAVCR2', 'PDCD1', 'LAG3', 'TIGIT', 'CTLA4', 'LAYN', 'ENTPD1'))
library(MySeuratWrappers)  
#需要展示的基因  
markers <- c('PRF1', 'IFNG', 'GNLY', 'NKG7', 'GZMB', 'GZMA', 'GZMH', 'KLRK1', 'KLRB1', 'KLRD1', 'CTSW', 'CST7')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(select_obj, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度


markers <- FindAllMarkers(select_obj, only.pos = TRUE)
write.table(markers, file = 'special_genes.csv', sep = ',')

saveRDS(select_obj, file = 'special_genes.rds')

#heatmap
allcolor <- c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080")
gene_features <- c('CXCL13', 'KRT86', 'TIGIT', 'CXCR6', 'CTLA4', 'SOX4', 'LAYN', 'TNFRSF18', 'DUSP4', 'PLPP1', 'IL7R', 'LMNA', 'GZMK', 'GPR183', 'ANXA1', 'XCL1', 'PIK3R1', 'MT2A', 'XCL2', 'HSPA6', 'SERPINH1', 'BAG3', 'HSPA1B', 'RHOB', 'HSPB1', 'DNAJB1', 'ZFAND2A', 'HSPA1A', 'STMN1', 'TYMS', 'UBE2C', 'RRM2', 'MKI67', 'KIAA0101', 'TOP2A', 'TK1', 'ASPM', 'TRAV1-2', 'IL4I1', 'ZBTB16', 'CA2', 'SLC4A10', 'IL23R', 'NCR3', 'RORC')
DoHeatmap(select_obj, features = gene_features, group.colors = allcolor)+scale_fill_gradientn(colors = c("#1A5B86","white","#921422"))#设置热图颜色  

FeaturePlot(select_obj, features = c('CD69', 'ITGAE', 'ITGA1', 'PDCD1', 'CTLA4', 'RUNX3', 'MKI67', 'CD38', 'ENTPD1'), min.cutoff = "q9")
FeaturePlot(select_obj, features = c('GZMB', 'GNLY', 'IFNG', 'CCL5', 'GZMA', 'CD3D', 'CD8A', 'NKG7', 'KLRG1'), min.cutoff = "q9")
FeaturePlot(select_obj, features = c('CD3E', "CXCR6", "ZNF683"), min.cutoff = "q3")

features_escc <- c("CD69", "ITGAE", "ITGA1", "CXCR6", "ZNF683","PDCD1", "CTLA4", "ENTPD1", "MKI67")
DotPlot(select_obj, features = features_escc) + RotatedAxis()

specolour=c("#FF3333", "#FFC5B2", '#34BB4F', "#336CB1", "#D3FF33")

new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4")
names(new.cluster.ids) <- levels(select_obj)
select_obj <- RenameIdents(select_obj, new.cluster.ids)
DimPlot(select_obj, reduction = "umap", label = TRUE, pt.size = 0.01, label.size = 6, raster=FALSE, cols = specolour)

#cd4 label clusters
c0 = c(0)
c1 = c(1)
c2 = c(2)
c3 = c(3)
c4 = c(4)

current.cluster.ids <- c(c0, c1, c2, c3, c4)

new.cluster.ids <- c(rep("C0",length(c1)), rep("C1",length(c1)), rep("C2",length(c2)), rep("C3",length(c3)), rep("C4",length(c4)))

select_obj@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(select_obj@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

saveRDS(select_obj, file = 'special_genes_labeled.rds')