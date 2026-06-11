setwd("E:\\论文\\拟南芥单细胞转录组学\\average_expression_correlation")
rds <- readRDS("E:\\论文\\拟南芥单细胞转录组学\\average_expression_correlation\\Early_flower_rename.slim.RDS")
proname <- "early_flower"
cluster <- "seurat_cluster"
library(Seurat)
library(ggplot2)
library(dplyr)
#统计每个cluster的细胞数量结果
cell.num <- data.frame(table(rds@meta.data$seurat_clusters))
write.table(cell.num,paste0("$proname",".cell.num.txt"))
#计算平均表达量
ave.exp <- AverageExpression(rds,group.by="seurat_clusters",layer="data")[["RNA"]]
colnames(ave.exp) <- as.vector(colnames(ave.exp))
write.table(ave.exp,paste0("$proname",".average.expression.txt"))
#读取注释信息
cell.ann <- read.csv("E:\\论文\\拟南芥单细胞转录组学\\average_expression_correlation\\earlyflower_ann.csv")

y <- as.character(cell.ann[,"Cluster"])
y <- paste0("g",y) 
ave.exp <- ave.exp[,y] 
colnames(ave.exp) <- paste(cell.ann[,"CellType"],cell.ann[,"Cluster"],sep="_")
ave.exp.dense <- as.matrix(ave.exp)
cor.exp <- as.data.frame(cor(ave.exp.dense))
cor.exp[["x"]] <- rownames(cor.exp)
cor.df <- tidyr::gather(data=cor.exp,y,correlation,-x) 
p1 <- ggplot(cor.df,aes(x,y,fill=correlation)) +  
  geom_point() + scale_fill_gradient(high="red",low = "white")
p1 <- p1 + theme(axis.text.x=element_text(angle=90,hjust = 1))
p1$"data"$"x" <- forcats::fct_inorder(p1$"data"$"x")
p1$"data"$"y" <- forcats::fct_inorder(p1$"data"$"y")
pdf(paste0("$proname", ".seurat.cluster.correlation.heatmap.pdf"), width = 9.55, height = 8)
p1
dev.off()

