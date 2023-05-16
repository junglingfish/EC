library(monocle)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)

setwd("D:/ZJU-FISH/BI/EC/data/GSE160269/")

cd010 <- readRDS("mye_des_labeled.rds")

cluster <- cd010$seurat_clusters
idx <- c(which(cluster==1), which(cluster==3), which(cluster==6), which(cluster==8))
cd010 <- cd010[,idx]
# expr_matrix <- as(as.matrix(cd010@assays$RNA@counts), 'sparseMatrix')
##��ȡ������Ϣ--ϸ����Ϣ(��������ϸ���ľ������ϸ�����ͼ�����Ϣ��ʵ����������Ϣ)
expr_matrix <- as(as.matrix(cd010@assays$RNA@counts), 'sparseMatrix')
##��ȡ������Ϣ��p_data(phenotype_data)���� 
p_data <- cd010@meta.data 
p_data$celltype <- cd010@active.ident  ##����ÿ��ϸ����ϸ��������Ϣ��p_data���档����Ѿ������򲻱��ظ�����
##��ȡ������Ϣ ���������͡�gc������
f_data <- data.frame(gene_short_name = row.names(cd010),row.names = row.names(cd010))
##expr_matrix��������f_data��������ͬ(gene number), expr_matrix��������p_data��������ͬ(cell number)

#����CDS����
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#��p_data��f_data��data.frameת��AnnotatedDataFrame����
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#����
cds <- detectGenes(cds, min_expr = 0.1) #��һ��������fData(cds)������һ��num_cells_expressed
print(head(fData(cds)))#��ʱ��13714������
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10)) 

##ʹ��seuratѡ��ĸ߱����???
express_genes <- VariableFeatures(cd010)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
#
# #ʹ��clusters����������
# deg.cluster <- FindAllMarkers(cd010)
# express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes)
# plot_ordering_genes(cds)
# 
# ##ʹ��monocleѡ��ĸ߱����???
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# cds <- setOrderingFilter(cds, disp.genes)
# plot_ordering_genes(cds)

# #��������
# diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1)
# #~�����Ǳ�ʾ��˭����������ı����������Ͽ���Ϊp_data����������
# head(diff)
# 
# ##������������Ϊ�켣�����Ļ���,��������ѡ���׼��qval<0.01,decreasing=F��ʾ����ֵ��������
# deg <- subset(diff, qval < 0.01) #ѡ��2829������
# deg <- deg[order(deg$qval,decreasing=F),]
# head(deg)
# 
# ##�������Ľ���ļ�����
# write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
# 
# ## �켣����������ӻ�
# ordergene <- rownames(deg)
# cds <- setOrderingFilter(cds, ordergene)
# #��һ���Ǻ���Ҫ�ģ������ǵõ���Ҫ�Ļ����б���������Ҫʹ��setOrderingFilter����Ƕ��cds���󣬺�����һϵ�в�����Ҫ���������list��
# #setOrderingFilter֮����Щ���򱻴�����cds@featureData@data[["use_for_ordering"]]������ͨ��table(cds@featureData@data[["use_for_ordering"]])�鿴
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
#??????ʹ��root_state��������������ʱ����ĸ������������ʱ����ɫͼ�п��Կ���������Ǹ�������stateͼ���Կ���������State1����Ҫ�����һ����Ϊ�������԰����²���
#cds <- orderCells(cds, root_state = 5) #��State5�����ʱ�������ʼ��

pdf("mye.pseudotime.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_tree = TRUE, show_backbone=TRUE) 
dev.off()

pdf("mye.state.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by="State", size=1, show_tree = TRUE, show_backbone=TRUE) 
dev.off()

pdf("mye.celltype.pdf",width = 10,height = 7)
plot_cell_trajectory(cds, color_by=Idents(cd010), size=1, show_tree = TRUE, show_backbone=TRUE) 
dev.off()