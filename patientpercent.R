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

mysce <- readRDS('special_genes.rds')

name <- Idents(mysce)
write.table(name, file = 'name_cluster_spe.xlsx', sep = ',')

# x <- read.csv("name_cluster_cd8.csv", sep = ',')
x <- read_excel("name_cluster_spe.xlsx", sheet = 1, col_names = FALSE)
nameList <- x$...1
nameList <- nameList[-1]

for(i in 1:length(nameList)){
  id1 = gregexpr("P", nameList[i])[[1]][1]
  id2 = gregexpr("I", nameList[i])[[1]][1]
  nameList[i] = substr(nameList[i], (id1+1), (id2-3))
}

a <- array(1:130, 130)
b <- array(0, 130)

for(i in 1:length(a)){
  count = 0
  flag = FALSE
  for(j in 1:length(nameList)){
    # chat = nameList[j]
    # print(i)
    # print(j)
    if(a[i] == as.numeric(nameList[j])){
      count = count + 1
      flag = TRUE
    }
  }
  if(flag == FALSE)
    b[i] = -1
  else
    b[i] = count/length(nameList)
}