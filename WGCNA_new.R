library(Seurat)
library(cowplot)
library(limma)
library(Matrix)
library(dplyr)
library('scran')
library(readxl)
library(monocle)
library(WGCNA)
library(OGFSC)
library(xlsx)
library(umap)
library(ggrepel)

setwd("D:/ZJU-FISH/BI/EC/data/")
seuratobj.data = readRDS('cd8.rds')


#WGCNA
cluster <- seuratobj.data$seurat_clusters
idx <- c(which(cluster==0), which(cluster==1), which(cluster==2), which(cluster==3), which(cluster==4), which(cluster==5), which(cluster==6))
seuratobj.data <- seuratobj.data[,idx]
#select cells
x <- 1:length(colnames(seuratobj.data@assays$RNA@counts))
id <- sample(x, size = 130, replace = F)

datExpr <- as.matrix(seuratobj.data@assays$RNA@counts)[,id]

#OGFSC
log2Data <- log2(datExpr +1)
## gene filtering by OGFSC
OGF <- OGFSC(log2Data, plot_option = 1, nBins = 30, minBinSize=300, LR_p=0.01,
             alpha=c(0.5), TW_threshold=0.0001) 
OGFSC_idx <- OGF$OGFSC_idx 
datExpr <- as.matrix(datExpr)[OGFSC_idx,]
datExpr <- t(datExpr)
datTraits <- as.matrix(seuratobj.data$seurat_clusters)
datTraits <- cbind(datTraits, seuratobj.data$orig.ident)
datTraits <- as.data.frame(datTraits[id,])
colnames(datTraits) <- c("cluster","ident")
for (i in 1:length(id)){
  datTraits$cluster[i] <- paste0("cluster", datTraits$cluster[i])
}
datTraits$cluster <- as.factor(datTraits$cluster)

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

# confirm beta
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate
par(mfrow = c(1,2))
cex1 = 0.9

# #power
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(datExpr, powerVector=powers,
#                         networkType="signed", verbose=5)
# 
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

net = blockwiseModules(datExpr, power = 1,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, maxBlockSize = 5000,
                       saveTOMFileBase = "300genes",
                       verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
pdf('Fig5d.pdf', width=6, height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# save
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

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
diag(dissTOM) = NA
# Call the plot function
table(moduleColors)
# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
        main = "Network heatmap plot, all genes")

#dataheatmap
nSamples <- nrow(datExpr)
nGenes = ncol(datExpr)
design=model.matrix(~0+ datTraits$cluster)
colnames(design)=levels(datTraits$cluster)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# sizeGrWindow(6,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
sta <- array(data = NA)
for (i in 1:nrow(moduleTraitCor)){
  sta[i] <- max(moduleTraitCor[i,])-min(moduleTraitCor[i,])
}
moduleTraitCor <- cbind(moduleTraitCor, sta)
moduleTraitCor <- moduleTraitCor[order(moduleTraitCor[,8], decreasing = T),]
# id <- which(moduleTraitCor[,7] < -0.13 | moduleTraitCor[,7] > 0.13)
# id <- id[-2]
# moduleTraitCor <- moduleTraitCor[id,-7]
moduleTraitCor <- moduleTraitCor[,-8]
textMatrix <- cbind(textMatrix, sta)
textMatrix <- textMatrix[order(textMatrix[,8], decreasing = T),]
textMatrix <- textMatrix[,-8]
# textMatrix <- textMatrix[id,]
# Display the correlation values within a heatmap plot
pdf('Fig5e.pdf', width=6, height=6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = T,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab.y = 0.3,
               cex.lab.x = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#gene
modNames = substring(colnames(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = 'p'))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
cluster1 = as.data.frame(design[,7])
names(cluster1) = "cluster1"
geneTraitSignificance = as.data.frame(cor(datExpr, cluster1, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(cluster1), sep="")
names(GSPvalue) = paste("p.GS.", names(cluster1), sep="")
module = "grey"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
module_gene <- abs(geneModuleMembership[moduleGenes, column])
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", colnames(cluster1)),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

x <- abs(geneModuleMembership[moduleGenes, column])
y <- abs(geneTraitSignificance[moduleGenes, 1])
data <- data.frame(x, y, row.names = rownames(geneModuleMembership)[moduleGenes])
data$name <- rownames(data)
data$name[which(data$x < 0.6 | data$y < 0.4)] <- NA
write.csv(data, "dot_data_c4_grey.csv")
data <- read.csv("dot_data.csv", row.names = 1)

tiff('Fig5f.tiff', units="in", width=9, height=8, res=300, compression = 'lzw')
ggplot(data)+
  geom_point(aes(x=x,y=y), size = 2, alpha=1, color="DarkCyan")+
  geom_text_repel(aes(x=x,y=y,label=name))+ xlab(paste("Module Membership in ", "grey", " module"))+
  ylab(paste("Gene significance for","cluster1"))+ theme_bw()+
  labs(title = paste("Module membership vs. gene significance\n", "cor=0.55,p=8.9e-24"))+ # need change
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("DarkCyan"))
dev.off()