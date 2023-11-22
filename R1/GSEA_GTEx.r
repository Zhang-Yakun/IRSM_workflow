setwd("E:/Rwork/lung/GTEx/")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 正常肺组织 -------------------------------------------------------------------

DEG<-read.csv("Lung(1).csv")

GSEA_df<-data.frame(SYMBOL=DEG$X, foldChange=DEG$logFC)

#a<-bitr(GSEA_df$ENTREZID,fromType = "ENTREZID",toType = "SYMBOL",OrgDb="org.Hs.eg.db")
#GSEA_df$ENTREZID<-GSEA_df$ENTREZID%>%as.character()
#GSEA_df2<-inner_join(a,GSEA_df)
#GSEA_df<-GSEA_df2[,-1]

sort.df<-GSEA_df[order(GSEA_df$foldChange,decreasing = TRUE),]
gene.expr<-sort.df$foldChange
names(gene.expr)<-sort.df$SYMBOL

head(gene.expr)

###GSEA
#BiocManager::install("clusterProfiler",force = TRUE)
library(AnnotationDbi)
library(GO.db)
library(clusterProfiler)
library(GSEABase)
library(enrichplot)

geneset<-clusterProfiler::read.gmt("E:/Rwork/lung/my_geneset.gmt")

egmt<-GSEA(gene.expr,TERM2GENE=geneset,pvalueCutoff = 1,verbose=FALSE)
head(egmt)

gseaplot2(egmt,geneSetID="my_geneset",pvalue_table = TRUE,subplots=c(1,2) )
write.csv(egmt,"GTEx_lung_gsea.csv",row.names=F)
