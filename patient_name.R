library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(tidyverse)
library(ROCR)
library(KernSmooth)
library(readxl)

setwd("D:/ZJU-FISH/BI/EC/data/")

mysce <- readRDS('special_genes_3clusters.rds')

name <- Idents(mysce)
write.table(name, file = 'name_cluster_spe_merged.csv', sep = ',')

# x <- read.csv("name_cluster_cd8.csv", sep = ',')
x <- read_excel("name_cluster_spe_merged.xlsx", sheet = 1, col_names = FALSE)
nameList <- x$...1
nameList <- nameList[-1]

for(i in 1:length(nameList)){
  id1 = gregexpr("P", nameList[i])[[1]][1]
  id2 = gregexpr("T", nameList[i])[[1]][1]
  nameList[i] = substr(nameList[i], (id1+1), (id2-1))
}

for(i in 1:length(nameList)){
  if(gregexpr("T", nameList[i])[[1]][1] > 0)
    nameList[i] = 1
  if(gregexpr("N", nameList[i])[[1]][1] > 0)
    nameList[i] = 0
}

TUM = c(1)
NOR = c(0)

current.sample.ids = c(TUM, NOR)
new.sample.ids = c(rep("Tumor",length(TUM)), rep("Normal",length(NOR)))

mysce@meta.data$sampletype <- plyr::mapvalues(x = as.integer(as.character(nameList)), from = current.sample.ids, to = new.sample.ids)

saveRDS(mysce, file = 'spe_sampletype_merged.rds')