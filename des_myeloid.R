library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(tidyverse)
library(ROCR)
library(KernSmooth)

setwd("D:/ZJU-FISH/BI/EC/data/GSE160269")

sce <- read.table(file="UMI_matrix_Myeloid.txt")

mysce <- CreateSeuratObject(counts = sce, project = "myeloid", min.cells = 3, min.features = 200)
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


mysce <- NormalizeData(mysce, normalization.method = "LogNormalize", scale.factor = 1e4) 
mysce <- FindVariableFeatures(mysce, selection.method = 'vst', nfeatures = 2000)
mysce <- ScaleData(mysce, vars.to.regress = "percent.mt")
mysce <- RunPCA(mysce, features = VariableFeatures(object = mysce))

# mysce <- JackStraw(mysce, num.replicate = 100)
# mysce <- ScoreJackStraw(mysce, dims = 1:20)
# ElbowPlot(mysce)

#10
mysce <- RunUMAP(mysce, reduction = "pca", dims = 1:10)
# mysce <- RunTSNE(mysce, reduction = "pca", dims = 1:15)
mysce <- FindNeighbors(mysce, dims = 1:10)
mysce <- FindClusters(mysce, resolution = 0.6)

p1 <- DimPlot(mysce, reduction = "umap")
p2 <- DimPlot(mysce, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#MC
features_escc <- c("CD163", 'CD68')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#M1-4
features_escc <- c("CD163", 'CD68', 'CD80', 'CD86', 'IRF5', 'STAT1', 'STAT2', 'STAT3', 'STAT6', 'IRF4', 'TGFB1', 'IL10', 'TNF', 'TLR2', 'ILR1')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#M2-9
features_escc <- c("CD163", 'CD68', 'MSR1', 'MRC1', 'CSF1R')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#MONO
features_escc <- c("FCN1", 'APOBEC3A', 'THBS1', 'CD14', 'FCGR3A', 'VCAN', 'FCGR2A', 'CSF1R')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#CLASSIC MONO-3
features_escc <- c("FCN1", 'APOBEC3A', 'THBS1', 'CD14', 'FCGR3A', 'VCAN', 'FCGR2A', 'CSF1R', 'CEBPD', 'CXCL14', 'MS4A6A', 'CLEC12A', 'S100A12', "FAM65B")
DotPlot(mysce, features = features_escc) + RotatedAxis()

#MEDIUM MONO-1
features_escc <- c("FCN1", 'APOBEC3A', 'THBS1', 'CD14', 'FCGR3A', 'VCAN', 'FCGR2A', 'CSF1R', 'ZFP36L2', 'CD300E', 'FAM65B')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#NOT CLASSIC MONO-0+5
features_escc <- c("FCN1", 'APOBEC3A', 'THBS1', 'CD14', 'FCGR3A', 'VCAN', 'FCGR2A', 'CSF1R', 'TCF7L2', 'KLF3', 'IKZF1', 'FLI1', 'ICAM2', 'CDH23', 'FCGR3B', 'CBL', 'CX3CR1', 'IFITM1', 'MTSS1', 'CDKN1C', 'SLC44A2', 'CD300E')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#DC
features_escc <- c("HLA-DQB2", 'HLA-DPB1', 'BIRC3', 'CD1C', 'FCER1A', 'CLEC4C')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#PDC-7
features_escc <- c("HLA-DQB2", 'HLA-DPB1', 'BIRC3', 'FCER1A', 'CLEC4C', 'SOX4', 'IRF4', 'IRF7', 'IRF8', 'LILRA4')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#CDC2-6
features_escc <- c("HLA-DQB2", 'HLA-DPB1', 'BIRC3', 'FCER1A', 'CLEC4C', 'CD1C', 'HMGA1', 'PFDA1', 'IRF4', 'ITGAX')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#CDC1-8
features_escc <- c("HLA-DQB2", 'HLA-DPB1', 'BIRC3', 'FCER1A', 'CLEC4C', 'BATF3', 'THBD', 'ID2', 'ETB3', 'CLEC9A', 'XCR1', 'CADM1')
DotPlot(mysce, features = features_escc) + RotatedAxis()

#MC2-2
features <- c("CPA3", "TPSB2", "TPSAB1")
DotPlot(mysce, features = features) + RotatedAxis()

# new.cluster.ids <- c("non-classical monocyte", "Intermediate monocyte", "Mast cell", "classical monocyte", "M1 macrophage", "non-classical monocyte", "cDC2", "pDC", "cDC1", "M2 macrophage")
# names(new.cluster.ids) <- levels(mysce)
# mysce <- RenameIdents(mysce, new.cluster.ids)
# allcolour=c("#BC3A23","#318937","#6F597B","#BD3B24","GREEN","#C28459", "#C1961D", "#B8BE77", "#C25D37")
# DimPlot(mysce, reduction = "umap", cols = allcolour, label = TRUE, pt.size = 0.5, label.size = 5)

new.cluster.ids <- c("M0", "M1", "M2", "M3", "M4", "M0", "M5", "M6", "M7", "M8")
names(new.cluster.ids) <- levels(mysce)
mysce <- RenameIdents(mysce, new.cluster.ids)
allcolour=c("#BC3A23","#318937","#6F597B","#BD3B24","GREEN","#C28459", "#C1961D", "#B8BE77", "#C25D37")
DimPlot(mysce, reduction = "umap", cols = allcolour, label = TRUE, pt.size = 0.5, label.size = 5)

#zhanbi
table(Idents(mysce))
prop.table(table(Idents(mysce)))

# library(ComplexHeatmap)
#heatmap
allcolor <- c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080")
gene_features <- c('HSPA6', 'SDS', 'BAG3', 'OLR1', 'HSPB1', 'ISG15', 'G0S2', 'FAM26F', 'PLIN2', 'ZFAND2A', 'CD14', 'FCGR3A', 'SPP1', 'FCGR2A', 'CD163', 'CCL18', 'TGFBI', 'MSR1', 'TPSAB1', 'TPSB2', 'CPA3', 'TPSD1', 'GATA2', 'CLU', 'MS4A2', 'LTC4S', 'NSMCE1', 'HDC', 'S100A8', 'S100A12', 'FCN1', 'NLRP3','APOBEC3A', 'VCAN', 'CD300E', 'LYZ', 'THBS1', 'CD14', 'CCL3', 'IL6', 'SPP1', 'C15orf48', 'IL1RN', 'IL1B', 'MARCO', 'EREG', 'HLA-DQB1', 'CD1E', 'CD1C', 'FCER1A', 'CLEC10A', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DPB1', 'IRF7', 'IRF4', 'LILRA4', 'PLD4', 'SOX4', 'IRF8', 'IL3RA', 'CCDC50', 'SPIB', 'CST7', 'IDO1', 'BIRC3', 'ANXA6', 'ID2', 'CD83', 'CD40', 'KLF4', 'KLF2', 'C1QC', 'C1QB', 'ATF3', 'MS4A4A', 'MRC1', 'CD163') 
DoHeatmap(mysce, features = gene_features, group.colors = allcolor)+scale_fill_gradientn(colors = c("#1A5B86","white","#921422"))#设置热图颜色  

# #violin
# VlnPlot(cd8, features = c('HAVCR2', 'PDCD1', 'LAG3', 'TIGIT', 'CTLA4', 'LAYN', 'ENTPD1'))
library(MySeuratWrappers)  
#需要展示的基因  
markers <- c('IL-33', 'TNF', 'TGFBI', 'IL15')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(mysce, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度


#result
features <- c("HLA-DPB1", "FCGR2A", "VCAN", "FCN1", "CD14", "FCGR3A", "ID2", "CD40", "CD1C", 'FCGR1A', "LILRA4", "CLEC4C", "CD68", "CD86", "TLR2", "CD163", "MRC1", 'MSR1', "CPA3", "TPSB2", "TPSAB1")
# DotPlot(mysce, features = features) + RotatedAxis()
DotPlot(mysce, features = features)+theme_bw()+theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066','#336699','#66CC66','#FFCC33'))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

#result
features <- c('CD69', 'ITGAE', 'ITGA1', 'CCR7', 'SELL', 'S1PR1')
# DotPlot(mysce, features = features) + RotatedAxis()
DotPlot(mysce, features = features)+theme_bw()+theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066','#336699','#66CC66','#FFCC33'))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 



mysce.markers <- FindAllMarkers(mysce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(mysce.markers, file = "mye_gene.csv", sep = ',')

# violin
# VlnPlot(cd8, features = c('HAVCR2', 'PDCD1', 'LAG3', 'TIGIT', 'CTLA4', 'LAYN', 'ENTPD1'))
library(MySeuratWrappers)  
#需要展示的基因  
markers <- c('IL10', 'TGFB1', 'C9orf26', 'IL15')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(mysce, features = markers, stacked=T,pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度

markers <- c('TNF', 'IL33', 'TGFB1', 'IL15')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(mysce, features = markers, stacked=T, pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度
#'CD274', 'PDCD1LG2', 'IDO1', 'HAVACR2', 'LGALS9', 'VISTA', 'NCR3LG1', 'HHLA2'

markers <- c('CD274', 'PDCD1LG2', 'IDO1', 'HAVCR2', 'LGALS9', 'VISTA', 'NCR3LG1', 'HHLA2')  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
VlnPlot(mysce, features = markers, stacked=T, pt.size=0, cols = my36colors, direction = "horizontal", x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())#不显示坐标刻度
#'CD274', 'PDCD1LG2', 'IDO1', 'HAVACR2', 'LGALS9', 'VISTA', 'NCR3LG1', 'HHLA2'

#label clusters
ncm = c(0, 5)
im = c(1)
mc = c(2)
cm = c(3)
m1 = c(4)
cdc2 = c(6)
pdc = c(7)
cdc1 = c(8)
m2 = c(9)

current.cluster.ids <- c(ncm, im, mc, cm, m1, cdc2, pdc, cdc1, m2)

new.cluster.ids <- c(rep("NCM",length(ncm)), rep("IM",length(im)), rep("MC",length(mc)), rep("CM",length(cm)), rep("M1",length(m1)), rep("CDC1",length(cdc1)), rep("PDC",length(pdc)), rep("CDC2",length(cdc2)), rep("M2",length(m2)))

mysce@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(mysce@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

saveRDS(mysce, file = 'mye_des_labeled.rds')

# #findmarker
# markers <- FindAllMarkers(mysce, only.pos = TRUE)
# write.table(markers, file = 'mye_markers.csv', sep = ',')



#mono/macro
mysce = readRDS('mye_des.rds')
mm = mysce[, Idents(mysce) %in% c( "mono/macro")]

mm <- NormalizeData(mm, normalization.method = "LogNormalize", scale.factor = 1e4) 
mm <- FindVariableFeatures(mm, selection.method = 'vst', nfeatures = 2000)
mm <- ScaleData(mm, vars.to.regress = "percent.mt")
mm <- RunPCA(mm, features = VariableFeatures(object = mm))

# mm <- JackStraw(mm, num.replicate = 100)
# mm <- ScoreJackStraw(mm, dims = 1:20)
# ElbowPlot(mm)

#11
mm <- RunUMAP(mm, reduction = "pca", dims = 1:13)
# mm <- RunTSNE(mm, reduction = "pca", dims = 1:15)
mm <- FindNeighbors(mm, dims = 1:13)
mm <- FindClusters(mm, resolution = 0.25)

# #12
# mm <- RunUMAP(mm, reduction = "pca", dims = 1:12)
# # mm <- RunTSNE(mm, reduction = "pca", dims = 1:15)
# mm <- FindNeighbors(mm, dims = 1:12)
# mm <- FindClusters(mm, resolution = 0.3)

p1 <- DimPlot(mm, reduction = "umap")
p2 <- DimPlot(mm, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#MC
features_escc <- c("CD163", 'CD68', 'ITGAX', 'CD1C', 'TGFB1', 'IL10', 'CD14', 'FCGR3A', 'CD80', 'CD86')
DotPlot(mm, features = features_escc) + RotatedAxis()

#MONO
features_escc <- c("FCN1", 'APOBEC3A', 'THBS1', 'CD14', 'FCGR3A', 'VCAN', 'FCGR2A', 'CSF1R', 'CEBPD', 'TCF7L2', 'KLF3', 'IKZF1', 'FLI1', 'ZFP36L2', 'CD300E')
DotPlot(mm, features = features_escc) + RotatedAxis()

#TAM
features_escc <- c("CCR2", 'CSF1R', 'MARCO', 'CD274', 'CD40', 'CCL2', 'CSF1', 'FCGR3A', 'PDGFB')
DotPlot(mm, features = features_escc) + RotatedAxis()

#dentric cell

new.cluster.ids <- c("MONO", "mono/macro", "mast cell", "dentric cell", "dentric cell", "dentric cell")
names(new.cluster.ids) <- levels(mysce)
mysce <- RenameIdents(mysce, new.cluster.ids)
DimPlot(mysce, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5, raster=FALSE)


