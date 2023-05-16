library(monocle)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)

setwd("D:/ZJU-FISH/BI/EC/data/")

cd010 <- readRDS("spe_nor_reset.rds")

cluster <- cd010$seurat_clusters
# idx <- c(which(cluster==1), which(cluster==3), which(cluster==6), which(cluster==8))
cd010 <- cd010[,idx]
# expr_matrix <- as(as.matrix(cd010@assays$RNA@counts), 'sparseMatrix')
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(cd010@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- cd010@meta.data 
p_data$celltype <- cd010@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(cd010),row.names = row.names(cd010))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#过滤
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10)) 

# ##使用seurat选择的高变基因???
# express_genes <- VariableFeatures(cd010)
# cds <- setOrderingFilter(cds, express_genes)
# plot_ordering_genes(cds)
#
# #使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(cd010)
# express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes)
# plot_ordering_genes(cds)
# 
##使用monocle选择的高变基因???
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

# #其他方法
# diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1)
# #~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
# head(diff)
# 
# ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
# deg <- subset(diff, qval < 0.01) #选出2829个基因
# deg <- deg[order(deg$qval,decreasing=F),]
# head(deg)
# 
# ##差异基因的结果文件保存
# write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
# 
# ## 轨迹构建基因可视化
# ordergene <- rownames(deg)
# cds <- setOrderingFilter(cds, ordergene)
# #这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
# #setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
# pdf("train.ordergenes.pdf")
# plot_ordering_genes(cds)
# dev.off()

# GM_state <- function(cds){
#   if (length(unique(pData(cds)$State)) > 1){
#     T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"0"]
#     return(as.numeric(names(T0_counts)[which
#                                        (T0_counts == max(T0_counts))]))
#   } else {
#     return (1)
#   }
# }


cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
# cds <- orderCells(cds, root_state = GM_state(cds))
#??????使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点

pdf("nor.pseudotime.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_tree = TRUE, show_backbone=TRUE)
dev.off()

pdf("nor.state.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by="State", size=1, show_tree = TRUE, show_backbone=TRUE) 
dev.off()

pdf("nor.celltype.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by=Idents(cd010), size=1, show_tree = TRUE, show_backbone=TRUE) 
dev.off()