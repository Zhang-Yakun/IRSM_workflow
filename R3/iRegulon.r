setwd("E:/Rwork/lung/TF/")

library(data.table)
library(dplyr)
library(stringr)
TF<-fread("TF.tsv")%>%as.data.frame()

##iRegulon的motif结果和track结果
motif_TF<-TF[1:78,]

track_TF<-TF[79:103,]
colnames(track_TF)<-track_TF[1,]
track_TF<-track_TF[-1,]

####读入表达数据计算组间差异

DEG<-read.csv("E:/Rwork/lung/limma_result_old.csv")

####track_TF
TF_deg<-subset(DEG,DEG$X%in%track_TF$`Transcription factor`)
TF_deg<-subset(TF_deg,TF_deg$P.Value<0.05)
TF_deg<-subset(TF_deg,TF_deg$logFC>0)
TF_deg$X

exp<-fread("exp_261.csv")%>%as.data.frame()
TF_exp<-exp[,TF_deg$X]
rownames(TF_exp)<-exp$V1

TF_exp$group<-exp$group

####motif_TF
motif1<-str_split(motif_TF$`Transcription factor`[1],",")%>%unlist
TF_deg1<-subset(DEG,DEG$X%in%motif1)

TF_deg1<-subset(TF_deg1,TF_deg1$P.Value<0.05)
TF_deg1$X

# 可视化差异转录因子 ---------------------------------------------------------------

library(ggstatsplot)

ggbetweenstats(data=TF_exp,x=group,y=MEF2A)
ggbetweenstats(data=TF_exp,x=group,y=NFKB1)
ggbetweenstats(data=TF_exp,x=group,y=MAX)
ggbetweenstats(data=TF_exp,x=group,y=MEF2C)
ggbetweenstats(data=TF_exp,x=group,y=STAT5A)


# 5个差异转录因子热图 --------------------------------------------------------------

TF_exp$sample<-rownames(TF_exp)
TF_exp$sum<-apply(TF_exp[,1:5],1,sum)
TF_exp2<-arrange(TF_exp,desc(group),sum)

group_list<-TF_exp2[,6:7]

exp_heat<-t(TF_exp2[,1:5])

identical(group_list$sample,colnames(exp_heat))

library(pheatmap)
cormat<-round(cor(exp_heat,method = "spearman"),2)

table(group_list$group)

annotation_col<-data.frame(
  Type = factor(group_list$group)) #按年龄分组
rownames(annotation_col) = colnames(exp_heat)

ann_colors=list(Type=c(young="#DDA678",old="#88BA9F"))

n <- t(scale(t(exp_heat)))#计算Z-score
n[n > 2] <- 2
n[n < -2] <- -2
df <- n

pheatmap(df,cellwidth = 1, cellheight = 16, fontsize = 14,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         #scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("#3178A3","white","#D86472"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors=ann_colors,
         treeheight_col = "0",#不画树
         border_color = "NA")

