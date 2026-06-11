library(DOSE)
library(org.At.tair.db)
library(topGO)
library(clusterProfiler)
library(gridExtra)
library(patchwork)
library(ggplot2)
setwd("E:/论文/拟南芥单细胞转录组学/GO_enrichmnet")
genes <- as.vector(read.csv("E:\\论文\\拟南芥单细胞转录组学\\GO_enrichmnet\\epi_deg_genes.csv",header=TRUE,sep="\t"))
prefix="S1_leaf"
gene_vector <- genes$geneID  
genelists=bitr(gene_vector,
                fromType="TAIR",
                toType="ENTREZID",
                OrgDb="org.At.tair.db")
id <- as.vector(genelists[,2])
ego.all <- enrichGO(
    gene=id,
    keyType ="ENTREZID",
    OrgDb=org.At.tair.db,
    ont="BP" ,
    pAdjustMethod="BH",
    pvalueCutoff=1,
    qvalueCutoff=1,
    readable=TRUE  
  )
pdf("GO Enrichment for s1_leaf_epi.pdf",width=7,height=5)
dotplot(ego.all,title=paste("GO Enrichment for",prefix,sep=""),orderBy="x")
dev.off()
