setwd("D:/workplace/mywork/")
# 表达谱预处理 ------------------------------------------------------------------

library(data.table)
library(dplyr)
expression <- fread("LUAD_readcount.genes.tpm.txt")%>%as.data.frame() #576个样本的基因表达谱

rownames(expression)<-expression$V1
expression<-expression[,-1]
expMatrix<-(log2(expression+1))

exp<-t(expMatrix)%>%as.data.frame()

sample_cli<-fread("LUAD_clinical_survival.csv")%>%as.data.frame()

cli_df<-subset(sample_cli,!is.na(sample_cli$age_at_diagnosis))
cli_df2<-data.frame(age=cli_df$age_at_diagnosis)
rownames(cli_df2)<-cli_df$sample
cli_df2$sample<-rownames(cli_df2)

intersect(rownames(exp),rownames(cli_df2))

exp$sample<-rownames(exp)

co_df<-inner_join(cli_df2,exp)

co_df$age%>%quantile()
# age 38~88

co_df$group<-ifelse(co_df$age>65,"old","young")
# old=289,young=268

rownames(co_df)<-co_df$sample

co_df<-co_df[,-1]

co_df2<-co_df%>%dplyr::select(group,everything())


# 差异分析 --------------------------------------------------------------------

group_list<- factor(co_df2$group,ordered = F) 

dta<-co_df2[,3:dim(co_df2)[2]] #87gene expression

limma_input<-t(dta)
limma_input<-as.data.frame(limma_input)

identical(colnames(limma_input),rownames(co_df2))


library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(limma_input , design)
contrast.matrix <- makeContrasts(old - young,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
DEG <- na.omit(tempOutput)

write.csv(DEG,"limma_result_old.csv")

# GSEA --------------------------------------------------------------------

DEG<-read.csv("limma_result_old.csv")

GSEA_df<-data.frame(SYMBOL=DEG$X, foldChange=DEG$logFC)

library(clusterProfiler)
library(org.Hs.eg.db)

df.id<-bitr(GSEA_df$SYMBOL,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

my.df<-merge(GSEA_df,df.id,by="SYMBOL",all=F)

sort.df<-my.df[order(my.df$foldChange,decreasing = TRUE),]
gene.expr<-sort.df$foldChange
names(gene.expr)<-sort.df$ENTREZID

head(gene.expr)

require(enrichplot)
require(clusterProfiler)

#######准备gmt文件-自己的基因集合


immune_gene90<-read.table("immunosense_gene.txt",header=F) #90个免疫衰老基因列表
immune_gene90<-immune_gene90$V1%>%unique()%>%as.character()

gset<-c("my_geneset","NA",immune_gene90)
gset<-gset%>%as.data.frame()%>%t()


write.table(gset,file="my_geneset.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

###GSEA
BiocManager::install("GSEA")
library(GSEA)
library(GSEABase)

geneset<-read.gmt(file.path(d,gmtfile))

egmt<-GSEA(gene.expr,TERM2GENE=geneset,verbose=FALSE)
head(egmt)