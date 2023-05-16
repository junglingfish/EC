library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ROCR)
library(KernSmooth)
library(immunarch)
library(readxl)

setwd("D:/ZJU-FISH/BI/EC/data/")

mysce <- readRDS('special_genes.rds')

# name <- Idents(mysce)
# write.table(name, file = 'name_cluster_spe.xlsx', sep = ',')

# x <- read.csv("name_cluster_cd8.csv", sep = ',')
x <- read_excel("name_cluster_spe.xlsx", sheet = 1, col_names = FALSE)
nameList <- x$...1
nameList <- nameList[-1]

# for(i in 1:length(nameList)){
#   id = gregexpr("I", nameList[i])[[1]][1]
#   nameList[i] = substr(nameList[i], (id+2), 50)
# }

all = array(126:130, 5)

for(i in 1:length(all)){
  seq = as.character(all[i])
  filepath = paste0("D:/ZJU-FISH/BI/EC/data/data_ZJU/data/filtered_contig_annotations.P", seq, "-N.csv")
  fileflag = file.exists(filepath)
  delete = array(0, 20000)
  count = 1
  marker = 1
  if(fileflag){
    # print(1)
    data <- read.table(filepath, sep=",", header=TRUE)
    buffer_barcode <- data[["barcode"]]
    for(j in 1:length(buffer_barcode)){
      buffer_barcode[j] = substr(buffer_barcode[j], 0, 16)
    }
    for(m in 1:length(buffer_barcode)){
      c = buffer_barcode[m]
      flag = FALSE
      for(n in 1:length(nameList)){
        id_P = gregexpr("P", nameList[n])[[1]][1]
        id_T = gregexpr("N", nameList[n])[[1]][1]
        patient_id = substr(nameList[n], (id_P+1), (id_T-1))
        if(patient_id == seq){id_P = gregexpr("P", nameList[n])[[1]][1]
            id_I = gregexpr("I", nameList[n])[[1]][1]
            patientname = substr(nameList[n], (id_I+2), (nchar(nameList[n])+1))
            if(patientname == c)
              flag = TRUE
        }
      }
      if(flag == FALSE){
        delete[count] = m
        count = count+1
      }
    }
    for(p in 1:length(delete)){
      if(delete[p] == 0){
        marker = p
        break
      }
    }
    a <- array(0, (marker-1))
    for(q in 1:(marker-1)){
      a[q] = delete[q]
    }
    data <- data[-a, ]
    write.csv(data, file = paste0("D:/ZJU-FISH/BI/EC/data/data_ZJU/data_spe/P", seq, "-N.csv"), sep = ',')
  }
}


### pic
setwd("D:/ZJU-FISH/BI/EC/data/data_ZJU/data_spe")

# immunarch
# immdata <- repLoad("filtered_contig_annotations.P1.csv")
immdata <- repLoad("D:/ZJU-FISH/BI/EC/data/data_ZJU/data_spe")

imm_rare <- repClonality(immdata$data, .method = "rare", .bound = c(1, 2))
imm_rare

vis(imm_rare)

write.table(imm_rare, file = 'immunarch_spe.csv', sep = ',')
# repExplore(immdata$data, "lens") %>% vis()
# 
# repClonality(immdata$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
# 
# # 可视化第一个克隆型的V基因分布图
# geneUsage(immdata$data[[1]]) %>% vis()  # Visualise the V-gene distribution for the first repertoire
# 
# # 构建不同克隆型之间共享的重叠克隆型热图
# repOverlap(immdata$data) %>% vis()  # Build the heatmap of public clonotypes shared between repertoires
# 
# # # 可视化克隆型的多样性
# # repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta)  # Visualise the Chao1 diversity of repertoires, grouped by the patient status

# immdata <- repLoad("D:/ZJU-FISH/BI/EC/data/data_ZJU/NT")
# 
# imm_rare <- repClonality(immdata$data, .method = "rare", .bound = c(1, 2))
# imm_rare

