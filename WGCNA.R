library(WGCNA)
library(dplyr)
library(Seurat)
library(tidyverse)
library(reshape2)
library(stringr)

setwd("D:/ZJU-FISH/BI/EC/data/")
cd8 = readRDS('mye_des_labeled.rds')

datadf <- as.matrix(cd8@assays$RNA@data)
idd1 <- cd8@meta.data
Inter.id1<-cbind(rownames(idd1),idd1$seurat_clusters)
rownames(Inter.id1)<-rownames(idd1)
colnames(Inter.id1)<-c("CellID","Celltype")
Inter.id1<-as.data.frame(Inter.id1)
head(Inter.id1)
Inter1<-datadf[,Inter.id1$CellID]
Inter2<-as.matrix(Inter1)
Inter2[1:4,1:4]

pseudocell.size = 100 ## 10 test
new_ids_list1 = list()
length(levels(factor(Inter.id1$Celltype)))

for (i in 1:length(levels(factor(Inter.id1$Celltype)))) {
  cluster_id = levels(factor(Inter.id1$Celltype))[i]
  cluster_cells <- rownames(Inter.id1[Inter.id1$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)     
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)    
  new_ids_list1[[i]] <- pseudo_ids      
}

new_ids <- unlist(new_ids_list1)
new_ids <- as.data.frame(new_ids)
head(new_ids)
new_ids_length <- table(new_ids)
new_ids_length

new_colnames <- rownames(new_ids)  ###add
#rm(all.data1)
gc()
colnames(datadf)  
all.data<-datadf[,as.character(new_colnames)] ###add
all.data <- t(all.data)###add
new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                    list(name=new_ids[,1]),FUN=mean)
rownames(new.data)<-new.data$name
new.data<-new.data[,-1]
new_ids_length<-as.matrix(new_ids_length)##

short<-which(new_ids_length<100)##
new_good_ids<-as.matrix(new_ids_length[-short,])##
result<-t(new.data)[,rownames(new_good_ids)]
dim(result)

# datExpr = GetAssayData(cd8, slot = "counts")
# write.table(datExpr, file = 'dataexpr.csv', sep = ',')
# 
# dataExpr = read.csv(file = 'dataexpr.csv')
# dataExpr <- as.data.frame(t(datExpr))

##gene
cd8 <- FindVariableFeatures(cd8, nfeatures = 800)
colnames(result)[grepl("[12]_Cel",colnames(result))]
Cluster1 <- result[intersect(Seurat::VariableFeatures(cd8),rownames(result)),]

###WGCNA
type = "unsigned"  # 官方推荐 "signed" 或 "signed hybrid"
corType = "pearson" # 相关性计算  官方推荐 biweight mid-correlation & bicor  corType: pearson or bicor 
corFnc = ifelse(corType=="pearson", cor, bicor)
corFnc
maxPOutliers = ifelse(corType=="pearson",1,0.05) # 对二元变量，如样本性状信息计算相关性时， # 或基因表达严重依赖于疾病状态时，需设置下面参数
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- as.matrix(Cluster1)

#gene
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

#trans
dataExpr <- as.data.frame(t(dataExprVar))
dim(dataExpr)
head(dataExpr)[,1:8]

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
gsg$goodSamples

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType="signed", verbose=5)

# par(mfrow = c(1,2))
# cex1 = 0.9
# # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# # 网络越符合无标度特征 (non-scale)
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",
#      ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"))
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red")
# # 筛选标准。R-square=0.85
# abline(h=0.85,col="red")
# 
# # Soft threshold与平均连通性
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
#      cex=cex1, col="red")
# power = sft$powerEstimate
# softPower  = power
# softPower

cor <- WGCNA::cor

net = blockwiseModules(dataExpr, power = 9, maxBlockSize = nGenes,#nGenes
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("dataExpr", ".tom"),
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleColors
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
head(MEs_col)
# ?plotEigengeneNetworks
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = T,
                      xLabelsAngle = 90)

###计算每个模块的特征向量基因，为某一特定模块第一主成分基因E。代表了该模块内基因表达的整体水平
dynamicColors <- c("#FF0000","#14B6FF","#32CD32","#FFA500","#FF6347","#DEB887", "#F08080")
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# 计算根据模块特征向量基因计算模块相异度：
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90) 

#heatmap
# which.module="turquoise"; 
# ME=MEs[, paste("ME",which.module, sep="")]
# par(mfrow=c(2,1), mar=c(0,4.1,4,2.05))
# plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),
#         nrgcols=30,rlabels=F,rcols=which.module,
#         main=which.module, cex.main=2)
# par(mar=c(2,2.3,0.5,0.8))
# barplot(ME, col=which.module, main="", cex.main=2,
#         ylab="eigengene expression",xlab="array sample")


load(net$TOMFiles, verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)
TOM[1:4,1:4]
#dim(TOM2)

dissTOM = 1-TOM
# dissTOM = TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
table(moduleColors)
# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
        main = "Network heatmap plot, all genes")

# TOMplot(plotTOM, net$dendrograms, moduleColors, 
#         main = "Network heatmap plot, all genes")

#later
newcd8<- CreateSeuratObject(Cluster1)

newcd8[["percent.mt"]] <- PercentageFeatureSet(newcd8, pattern = "^CD")

newcd8 <- FindVariableFeatures(newcd8, selection.method = "vst", nfeatures = 2000)
newcd8 %>% NormalizeData( normalization.method = "LogNormalize", scale.factor = 10000)%>%
  FindVariableFeatures( selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features=VariableFeatures(newcd8),vars.to.regress = "percent.mt")  %>%
  RunPCA(features = VariableFeatures(object = newcd8)) %>%
  FindNeighbors( dims = 1:10) %>%
  FindClusters( resolution = 0.5) %>%
  BuildClusterTree() %>%
  RunUMAP( dims = 1:10)  -> newcd8

head(newcd8@meta.data)


moduleTraitCor_noFP <- cor(mergedMEs, newcd8@meta.data, use = "p");
moduleTraitPvalue_noFP = corPvalueStudent(moduleTraitCor_noFP, nSamples); 
textMatrix_noFP <- paste(signif(moduleTraitCor_noFP, 2), "\n(", signif(moduleTraitPvalue_noFP, 1), ")", sep = ""); 
par(mar = c(10, 8.5, 3, 3)); 
labeledHeatmap(Matrix = moduleTraitCor_noFP, 
               xLabels = names(newcd8@meta.data), 
               yLabels = names(mergedMEs), 
               ySymbols = names(mergedMEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix_noFP,
               setStdMargins = FALSE, 
               cex.text = 0.65, 
               zlim = c(-1,1), 
               main = paste("Module-trait relationships")) 