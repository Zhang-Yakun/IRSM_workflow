
# GSE1 GSEA--------------------------------------------------------------------

setwd("E:/Rwork/lung/GSE40419/")

DEG<-read.csv("nrDEG_GSE40419.csv")

GSEA_df<-data.frame(SYMBOL=DEG$X, foldChange=DEG$logFC)

library(clusterProfiler)
library(org.Hs.eg.db)

df.id<-bitr(GSEA_df$SYMBOL,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

my.df<-merge(GSEA_df,df.id,by="SYMBOL",all=F)

sort.df<-my.df[order(my.df$foldChange,decreasing = TRUE),]
gene.expr<-sort.df$foldChange
names(gene.expr)<-sort.df$ENTREZID

head(gene.expr)

library(enrichplot)

###GSEA
#BiocManager::install("clusterProfiler",force = TRUE)
library(AnnotationDbi)
library(GO.db)
library(clusterProfiler)
library(GSEABase)
library(enrichplot)

geneList=GSEA_df$foldChange
names(geneList)=GSEA_df$SYMBOL
geneList=sort(geneList,decreasing = T)

geneset<-clusterProfiler::read.gmt("E:/Rwork/lung/my_geneset.gmt")

egmt<-GSEA(geneList,TERM2GENE=geneset,verbose=FALSE)
head(egmt)

gseaplot2(egmt,geneSetID="my_geneset",pvalue_table = TRUE,subplots=c(1,2) )
write.csv(egmt,"GSEA_GSE40419.csv",row.names=F)

# GSE2 GSEA 富集结果不好--------------------------------------------------------------------

setwd("E:/Rwork/lung/GSE81089/")

DEG<-read.csv("nrDEG_GSE81089.csv")

GSEA_df<-data.frame(SYMBOL=DEG$X, foldChange=DEG$logFC)

library(clusterProfiler)
library(org.Hs.eg.db)

df.id<-bitr(GSEA_df$SYMBOL,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

my.df<-merge(GSEA_df,df.id,by="SYMBOL",all=F)

sort.df<-my.df[order(my.df$foldChange,decreasing = TRUE),]
gene.expr<-sort.df$foldChange
names(gene.expr)<-sort.df$ENTREZID

head(gene.expr)

library(enrichplot)

###GSEA
#BiocManager::install("clusterProfiler",force = TRUE)
library(AnnotationDbi)
library(GO.db)
library(clusterProfiler)
library(GSEABase)
library(enrichplot)

geneList=GSEA_df$foldChange
names(geneList)=GSEA_df$SYMBOL
geneList=sort(geneList,decreasing = T)

geneset<-clusterProfiler::read.gmt("E:/Rwork/lung/my_geneset.gmt")

egmt<-GSEA(geneList,TERM2GENE=geneset,verbose=FALSE)
head(egmt)

gseaplot2(egmt,geneSetID="my_geneset",pvalue_table = TRUE,subplots=c(1,2) )
write.csv(egmt,"GSEA_GSE81089.csv",row.names=F)

# GSE3 GSEA--------------------------------------------------------------------

####GSEA输入
setwd("E:/Rwork/lung/GSE209891/")

DEG3<-read.csv("nrDEG_GSE209891.csv")

GSEA_df3<-data.frame(SYMBOL=DEG3$GeneName, foldChange=DEG3$logFC)

library(clusterProfiler)
library(org.Mm.eg.db)

df.id3<-bitr(GSEA_df3$SYMBOL,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

my.df3<-merge(GSEA_df3,df.id3,by="SYMBOL",all=F)

sort.df3<-my.df3[order(my.df3$foldChange,decreasing = TRUE),]
gene.expr3<-sort.df3$foldChange
names(gene.expr3)<-sort.df3$ENTREZID

head(gene.expr3)

library(enrichplot)

###读入基因集
#BiocManager::install("clusterProfiler",force = TRUE)
library(AnnotationDbi)
library(GO.db)
library(clusterProfiler)
library(GSEABase)
library(enrichplot)

geneList3=GSEA_df3$foldChange
names(geneList3)=GSEA_df3$SYMBOL
geneList3=sort(geneList3,decreasing = T)

geneset<-clusterProfiler::read.gmt("E:/Rwork/lung/my_geneset.gmt")

######将鼠基因转换为人基因
library(biomaRt)

listDatasets(useMart("ensembl"))#查看ensembl数据库中可用的生物数据
human<-useMart("ensembl","hsapiens_gene_ensembl") #检索人类数据集
human1<-useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
mouse<-useMart("ensembl","mmusculus_gene_ensembl") #检索鼠数据集
mouse1<-useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "asia")
#有效的镜像选项是“ www”，“ uswest”，“ useast”，“ asia”。如果未指定镜像，则将使用www.ensembl.org上的主站点

listAttributes(human) #查看可用的参数
#mgi_symbol表示小鼠的基因名，hgnc_symbol表示人基因名

#转换ID,mouse2human
getLDS(
  values = names(geneList3),mart = mouse1,
  attributes = "mgi_symbol",filters = "mgi_symbol",
  martL = human1,
  attributesL = c("hgnc_symbol","chromosome_name"))
####直接改变基因大小写

head(geneset$gene)

library(stringr)
geneset$gene<-str_to_title(geneset$gene)

#library(Hmisc);capitalize(tolower(geneset$gene))

###GSEA
egmt3<-GSEA(geneList3,TERM2GENE=geneset,verbose=FALSE)
head(egmt3)

gseaplot2(egmt3,geneSetID="my_geneset",pvalue_table = TRUE,subplots=c(1,2) )
write.csv(egmt3,"GSEA_GSE209891.csv",row.names=F)
