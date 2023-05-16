library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ROCR)
library(KernSmooth)
library(immunarch)

setwd("D:/ZJU-FISH/BI/EC/data/data_ZJU/")

# immunarch
# immdata <- repLoad("filtered_contig_annotations.P1.csv")
immdata <- repLoad("D:/ZJU-FISH/BI/EC/data/data_ZJU/data")

imm_rare <- repClonality(immdata$data, .method = "rare", .bound = c(1, 2))
imm_rare

vis(imm_rare)

write.table(imm_rare, file = 'immunarch_all.csv', sep = ',')
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

immdata <- repLoad("D:/ZJU-FISH/BI/EC/data/data_ZJU/NT")

imm_rare <- repClonality(immdata$data, .method = "rare", .bound = c(1, 2))
imm_rare

vis(imm_rare)

write.table(imm_rare, file = 'immunarch_nt.csv', sep = ',')