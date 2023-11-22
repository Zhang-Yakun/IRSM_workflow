
setwd("E:/Rwork/lung/")
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

write.csv(co_df,"LUAD_tpm_age.csv",row.names = F)

#####剔除正常样本

library(stringr)

co_df$patient<-ifelse(co_df$sample%>%str_sub(14,15)=="11","normal","LUAD")
co_df<-co_df%>%dplyr::select(patient,everything())
table(co_df$patient)
# LUAD=498,normal=59

co_df$group<-ifelse(co_df$age<=59,"young",ifelse(co_df$age>=73,"old","middle"))

rownames(co_df)<-co_df$sample
co_df<-co_df[,-3]

co_df2<-subset(co_df,co_df$patient=="LUAD")

co_df2$age%>%quantile()
co_df2<-co_df2%>%dplyr::select(group,everything())

co_df3<-subset(co_df2,co_df2$group%in%c("young","old"))
##old=125,young=136

write.csv(co_df3,"exp_261.csv")

# 差异分析 --------------------------------------------------------------------

group_list<- factor(co_df3$group,ordered = F) 

dta<-co_df3[,4:dim(co_df3)[2]] #87gene expression

limma_input<-t(dta)
limma_input<-as.data.frame(limma_input)

identical(colnames(limma_input),rownames(co_df3))


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

DEG<-read.csv("E:/Rwork/lung/limma_result_old.csv")

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
BiocManager::install("clusterProfiler",force = TRUE)
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
