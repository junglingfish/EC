library(Seurat)
library(cowplot)
library(limma)
library(Matrix)
library(umap)
library(dplyr)
library('scran')
library(readxl)
library(pheatmap)
library(RColorBrewer) 
library(monocle)
library(ggpubr)
library(tidyverse)

#cd8 data
setwd("D:/ZJU-FISH/BI/EC/data/")

seuratobj.data <- readRDS("cd4.rds")

# Select the clusters want to map
cluster = seuratobj.data$seurat_clusters
idx_cluster <- c(which(cluster==0), which(cluster==1), which(cluster==2), which(cluster==3), which(cluster==4), which(cluster==5), which(cluster==6), which(cluster==7))
seuratobj.sel <- seuratobj.data[,idx_cluster]

# cell cycle, dotplot
DefaultAssay(seuratobj.data) <- "RNA"
seuratobj.sel <- seuratobj.data[,idx_cluster]

# score part
setwd("D:/ZJU-FISH/BI/EC/data/score/")

# # cyto score
# T <- read_excel("score_cent.xlsx", sheet = 1, col_names = FALSE)
# geneList <- T$...1
# geneList <- geneList[-1]
# idx_cyto <- match(geneList, rownames(seuratobj.data))
# idx_cyto <- idx_cyto[ !is.na(idx_cyto) ]

# exha list
T <- read_excel("score_cent.xlsx", sheet = 2, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_ex <- match(geneList, rownames(seuratobj.data))
idx_ex <- idx_ex[ !is.na(idx_ex) ]

#naive score
T <- read_excel("score_cent.xlsx", sheet = 3, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_na <- match(geneList, rownames(seuratobj.data))
idx_na <- idx_na[ !is.na(idx_na) ]

#treg score
T <- read_excel("score_cent.xlsx", sheet = 4, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_tr <- match(geneList, rownames(seuratobj.data))
idx_tr <- idx_tr[ !is.na(idx_tr) ]

#mine_trm score
T <- read_excel("score_mine.xlsx", sheet = 1, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
geneList <- toupper(geneList)
idx_mitrm <- match(geneList, rownames(seuratobj.data))
idx_mitrm <- idx_mitrm[ !is.na(idx_mitrm) ]

#all score
# expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_cyto,])
# cyto_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_ex,])
ex_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_na,])
na_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_tr,])
tr_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_mitrm,])
mitrm_scores <- colMeans(expression, na.rm = FALSE)

expr <- rbind(ex_scores, na_scores, tr_scores, mitrm_scores)
rownames(expr) <- c("exha", 'naive', 'treg', 'trm')

#calculate
sd <- colSds(t(expr))*sqrt((dim(expr)[2]-1)/(dim(expr)[2]))
mean <- colMeans(t(expr), na.rm = FALSE, dims = 1)
z_score <- matrix(0, dim(expr)[1], dim(expr)[2])
for (i in 1:4){
  z_score[i,] <- (expr[i,]-mean[i])/sd[i]
}
rownames(z_score) <- rownames(expr)
colnames(z_score) <- colnames(expr)
temp <- rbind(seuratobj.sel@assays$RNA@data, z_score)
seuratobj.sel@assays$RNA@data <- temp

markers.to.plot <- c("exha", 'naive', 'treg', 'trm')
tiff('cd4_score.tiff', units="in", width=8, height=6.5, res=300, compression = 'lzw')
DotPlot(seuratobj.sel, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 16) + RotatedAxis() +  ggtitle("CD4+ T cell cycle")
dev.off()


# #effect
# T <- read_excel("score_eff.xlsx", sheet = 1, col_names = FALSE)
# geneList <- T$...1
# geneList <- geneList[-1]
# geneList <- toupper(geneList)
# idx_eff <- match(geneList, rownames(seuratobj.data))
# idx_eff <- idx_eff[ !is.na(idx_eff) ]
# 
# #TRM score
# T <- read_excel("TRM-TCM.xlsx", sheet = 1, col_names = FALSE)
# geneList <- T$...1
# geneList <- geneList[-1]
# geneList <- toupper(geneList)
# idx_trm <- match(geneList, rownames(seuratobj.data))
# idx_trm <- idx_trm[ !is.na(idx_trm) ]
# 
# #mem score
# T <- read_excel("mem.xlsx", sheet = 1, col_names = FALSE)
# geneList <- T$...1
# geneList <- geneList[-1]
# geneList <- toupper(geneList)
# idx_mem <- match(geneList, rownames(seuratobj.data))
# idx_mem <- idx_mem[ !is.na(idx_mem) ]
# 
# #all score
# expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_eff,])
# eff_scores <- colMeans(expression, na.rm = FALSE)
# expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_trm,])
# trm_scores <- colMeans(expression, na.rm = FALSE)
# expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_tcm,])
# tcm_scores <- colMeans(expression, na.rm = FALSE)
# expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_mem,])
# mem_scores <- colMeans(expression, na.rm = FALSE)
# 
# expr <- rbind(eff_scores, trm_scores, tcm_scores, mem_scores)
# rownames(expr) <- c("eff", "trm", 'tcm', 'mem')
# 
# #calculate
# sd <- colSds(t(expr))*sqrt((dim(expr)[2]-1)/(dim(expr)[2]))
# mean <- colMeans(t(expr), na.rm = FALSE, dims = 1)
# z_score <- matrix(0, dim(expr)[1], dim(expr)[2])
# for (i in 1:4){
#   z_score[i,] <- (expr[i,]-mean[i])/sd[i]
# }
# rownames(z_score) <- rownames(expr)
# colnames(z_score) <- colnames(expr)
# temp <- rbind(seuratobj.sel@assays$RNA@data, z_score)
# seuratobj.sel@assays$RNA@data <- temp
# 
# markers.to.plot <- c("eff", "trm", 'tcm', 'mem')
# tiff('cd4_score_2.tiff', units="in", width=8, height=6.5, res=300, compression = 'lzw')
# DotPlot(seuratobj.sel, features = rev(markers.to.plot), cols = c("grey", "blue"), dot.scale = 16) + RotatedAxis() +  ggtitle("CD4+ T cell cycle")
# dev.off()

# cell cycle, dotplot
# G1 list
T <- read_excel("cell cycle genes.xlsx", sheet = 1, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_G1 <- match(geneList, rownames(seuratobj.data))
idx_G1 <- idx_G1[ !is.na(idx_G1) ]
# S list
T <- read_excel("cell cycle genes.xlsx", sheet = 2, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_S <- match(geneList, rownames(seuratobj.data))
idx_S <- idx_S[ !is.na(idx_S) ]
# G2 list
T <- read_excel("cell cycle genes.xlsx", sheet = 3, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_G2 <- match(geneList, rownames(seuratobj.data))
idx_G2 <- idx_G2[ !is.na(idx_G2) ]
# M list
T <- read_excel("cell cycle genes.xlsx", sheet = 4, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_M <- match(geneList, rownames(seuratobj.data))
idx_M <- idx_M[ !is.na(idx_M) ]

expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_G1,])
G1_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_S,])
S_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_G2,])
G2_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_M,])
M_scores <- colMeans(expression, na.rm = FALSE)
expr <- rbind(G1_scores, S_scores, G2_scores, M_scores)
rownames(expr) <- c("G1", "S", "G2", 'M')

sd <- colSds(t(expr))*sqrt((dim(expr)[2]-1)/(dim(expr)[2]))
mean <- colMeans(t(expr), na.rm = FALSE, dims = 1)
z_score <- matrix(0, dim(expr)[1], dim(expr)[2])
for (i in 1:4){
  z_score[i,] <- (expr[i,]-mean[i])/sd[i]
}
rownames(z_score) <- rownames(expr)
colnames(z_score) <- colnames(expr)
temp <- rbind(seuratobj.sel@assays$RNA@data, z_score)
seuratobj.sel@assays$RNA@data <- temp



markers.to.plot <- c("G1", "S", "G2", "M")
tiff('cd8_circle.tiff', units="in", width=6, height=4.5, res=300, compression = 'lzw')
DotPlot(seuratobj.sel, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 16) + RotatedAxis() +  ggtitle("NK+ cell cycle")
dev.off()
