library(Seurat)
library(readxl)
library(pheatmap)
library(RColorBrewer) 
library(ggplot2)
library(dplyr)
library(tidyverse)

rm(list=ls())
graphics.off()

setwd("D:/ZJU-FISH/BI/EC/data/")

cd8 <- readRDS(file = 'cd8.rds')

cd4 <- readRDS(file = 'cd4.rds')

#cd8
cluster8 <- cd8$seurat_clusters
idx_C0 <- which(cluster8==0)
cd8.C0 <- cd8[,idx_C0]
idx_C1 <- which(cluster8==1)
cd8.C1 <- cd8[,idx_C1]
idx_C2 <- which(cluster8==2)
cd8.C2 <- cd8[,idx_C2]
idx_C3 <- which(cluster8==3)
cd8.C3 <- cd8[,idx_C3]
idx_C4 <- which(cluster8==4)
cd8.C4 <- cd8[,idx_C4]
idx_C5 <- which(cluster8==5)
cd8.C5 <- cd8[,idx_C5]
idx_C6 <- which(cluster8==6)
cd8.C6 <- cd8[,idx_C6]

C8 <- c("C0", "C1","C2", "C3", "C4", "C5", "C6")
Cn8 <- c(0, 1, 2, 3, 4, 5, 6)

cd8.sub <- list()
for (i in 1:length(Cn8)){
  idx <- which(cluster8==Cn8[i])
  cd8.sub[[i]] <- cd8[,idx]
  #  assign(paste0("TCR_barcodes_C", Cn[i]), list())
  #  assign(paste0("TCR_cloneTypes_C", Cn[i]), list())
}

#RESTART
TCR_barcode_buffer = NULL
TCR_cloneType_buffer = NULL

setwd("D:/ZJU-FISH/BI/EC/data/data_ZJU/data")


### load  sampleA-tumor TCR results
data <- read.table("./filtered_contig_annotations.P130-T.csv", sep=",", header=TRUE)
buffer_barcode <- data[["barcode"]]
buffer_cloneType <- data[["raw_clonotype_id"]]
buffer_barcode_unique <- unique(buffer_barcode)
idx_selected <- match(buffer_barcode_unique, buffer_barcode)
buffer_barcode <- buffer_barcode[idx_selected]
buffer_cloneType <- buffer_cloneType[idx_selected]
idx_none <- which(buffer_cloneType=="None")
buffer_barcode <- buffer_barcode[-idx_none]
buffer_cloneType <- buffer_cloneType[-idx_none]
buffer_barcode <- sub("-.*", "", buffer_barcode) # extract the barcode sequence
buffer_barcode <- paste0('P130T.I.', buffer_barcode)
TCR_barcode_buffer = c(TCR_barcode_buffer, buffer_barcode)
TCR_cloneType_buffer = c(TCR_cloneType_buffer, buffer_cloneType)

for (i in 1:length(Cn8)){
  # idx1 <- which(cd8.sub[[i]]$sample == "tumor")
  # idx2 <- which(cd8.sub[[i]]$group == "sampleA")
  # idx_sample <- intersect(idx1, idx2)
  idx_sample <- which(cd8.sub[[i]]$orig.ident == "escc")
  Cell_barcodes_buffer <- colnames(cd8.sub[[i]])[idx_sample]
  Cell_barcodes_buffer <- sub("_.*", "", Cell_barcodes_buffer)
  idx_TCRinRNAseq <- match(Cell_barcodes_buffer, buffer_barcode)
  idx_NA <- which(is.na(idx_TCRinRNAseq))
  idx_TCRinRNAseq <- idx_TCRinRNAseq[-idx_NA]
  assign(paste0("TCR_barcodes_C8_C", Cn8[i], "_P130_T_tumor"), buffer_barcode[idx_TCRinRNAseq])
  assign(paste0("TCR_cloneTypes_C8_C", Cn8[i], "_P130_T_tumor"), buffer_cloneType[idx_TCRinRNAseq])
}

#cd4
cluster4 <- cd4$seurat_clusters
idx_C0 <- which(cluster4==0)
cd4.C0 <- cd4[,idx_C0]
idx_C1 <- which(cluster4==1)
cd4.C1 <- cd4[,idx_C1]
idx_C2 <- which(cluster4==2)
cd4.C2 <- cd4[,idx_C2]
idx_C3 <- which(cluster4==3)
cd4.C3 <- cd4[,idx_C3]
idx_C4 <- which(cluster4==4)
cd4.C4 <- cd4[,idx_C4]
idx_C5 <- which(cluster4==5)
cd4.C5 <- cd4[,idx_C5]
idx_C6 <- which(cluster4==6)
cd4.C6 <- cd4[,idx_C6]
idx_C7 <- which(cluster4==7)
cd4.C7 <- cd4[,idx_C7]

C4 <- c("C0", "C1","C2", "C3", "C4", "C5", "C6", "C7")
Cn4 <- c(0, 1, 2, 3, 4, 5, 6, 7)

cd4.sub <- list()
for (i in 1:length(Cn4)){
  idx <- which(cluster4==Cn4[i])
  cd4.sub[[i]] <- cd4[,idx]
  #  assign(paste0("TCR_barcodes_C", Cn[i]), list())
  #  assign(paste0("TCR_cloneTypes_C", Cn[i]), list())
}

for (i in 1:length(Cn4)){
  # idx1 <- which(cd8.sub[[i]]$sample == "tumor")
  # idx2 <- which(cd8.sub[[i]]$group == "sampleA")
  # idx_sample <- intersect(idx1, idx2)
  idx_sample <- which(cd4.sub[[i]]$orig.ident == "escc")
  Cell_barcodes_buffer <- colnames(cd4.sub[[i]])[idx_sample]
  Cell_barcodes_buffer <- sub("_.*", "", Cell_barcodes_buffer)
  idx_TCRinRNAseq <- match(Cell_barcodes_buffer, buffer_barcode)
  idx_NA <- which(is.na(idx_TCRinRNAseq))
  idx_TCRinRNAseq <- idx_TCRinRNAseq[-idx_NA]
  assign(paste0("TCR_barcodes_C4_C", Cn4[i], "_P130_T_tumor"), buffer_barcode[idx_TCRinRNAseq])
  assign(paste0("TCR_cloneTypes_C4_C", Cn4[i], "_P130_T_tumor"), buffer_cloneType[idx_TCRinRNAseq])
}


library(xlsx)
wb <- xlsx::createWorkbook()
# sheets <- getSheet()
sheet1 <- createSheet(wb,sheetName = "clone")
sheet2 <- createSheet(wb,sheetName = "cd8")
sheet3 <- createSheet(wb,sheetName = "cd4")
cols1 = 1
cols2 = 1
cols3 = 1


#pic cd8
library(gridExtra)
library(dplyr)
Cn8 <- c(0, 1, 2, 3, 4, 5, 6)
C8 <- c("C0", "C1","C2", "C3", "C4", "C5", "C6")
sample <- c("P130_T-tumor")
buffer8 <- matrix(0, length(sample), length(Cn8))
for (i in 1:length(Cn8)){
  temp <- eval(parse(text = paste0("TCR_cloneTypes_C8_C", Cn8[i], "_P130_T_tumor")))
  flag <- paste0('P130_T_CD8_C', Cn4[i])
  temp <- append(temp, flag)
  temps <- rev(temp)
  addDataFrame(temps,sheet1,row.names = F, col.names = F, startColumn = cols1)
  cols1 = cols1 + 1
  a <- table(temp)
  buffer8[1,i] <- sum(a[which(a>1)])/sum(a)
}

#pic cd4
library(gridExtra)
library(dplyr)
Cn4 <- c(0, 1, 2, 3, 4, 5, 6, 7)
C4 <- c("C0", "C1","C2", "C3", "C4", "C5", "C6", "C7")
sample <- c("P130_T-tumor")
buffer4 <- matrix(0, length(sample), length(Cn4))
for (i in 1:length(Cn4)){
  temp <- eval(parse(text = paste0("TCR_cloneTypes_C4_C", Cn4[i], "_P130_T_tumor")))
  flag <- paste0('P130_T_CD4_C', Cn4[i])
  temp <- append(temp, flag)
  temps <- rev(temp)
  addDataFrame(temps,sheet1,row.names = F, col.names = F, startColumn = cols1)
  cols1 = cols1 + 1
  a <- table(temp)
  buffer4[1,i] <- sum(a[which(a>1)])/sum(a)
}

buffer8 <- append(buffer8, 'P130_T_cd8')
addDataFrame(buffer8,sheet2,row.names = F, col.names = F, startColumn = cols2)
cols2 = cols2 + 1
buffer4 <- append(buffer4, 'P130_T_cd4')
addDataFrame(buffer4,sheet3,row.names = F, col.names = F, startColumn = cols3)
cols3 = cols3 + 1
# buffer8
# buffer4


setwd("D:/ZJU-FISH/BI/EC/data/data_ZJU/result/cloneresult")
saveWorkbook(wb,file = "clone_P130_T.xlsx")