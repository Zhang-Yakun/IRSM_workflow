
library(car)
library(sva)
library(genefilter)
library(preprocessCore)
library(ridge)
library(pRRophetic)

setwd("E:/Rwork/lung/R5/drug/")

library(ggplot2)
data(bortezomibData)
dim(exprDataBortezomib)
exprDataBortezomib[1:4,1:4]
boxplot(exprDataBortezomib[,1:4])
# 264个样品的表达量矩阵
table(studyResponse)

predictedPtype <- pRRopheticPredict(testMatrix=exprDataBortezomib, 
                                    drug="Bortezomib",
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2014")
library(ggplot2)
data(bortezomibData)
dim(exprDataBortezomib)
exprDataBortezomib[1:4,1:4]
boxplot(exprDataBortezomib[,1:4])
# 264个样品的表达量矩阵
table(studyResponse)

predictedPtype <- pRRopheticPredict(testMatrix=exprDataBortezomib, 
                                    drug="Bortezomib",
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2014")


# GDSC --------------------------------------------------------------------

setwd("E:/Rwork/lung/R5/drug/")

library(dplyr)
library(data.table)

d1<-fread("GDSC1_fitted_dose_response_25Feb20.csv")%>%as.data.frame()
##### GDSC2包含GDSC1，GDSC1是该网站上可用的原始数据集(2009-2015年间收集)的扩展。
#####而GDSC2则基于改进的技术、设备和程序等所得的最新的数据(2015-至今)。例如：GDSC1使用DNA染料(Syto60)，而GDSC2使用代谢测定法(Resazurin / CellTiter-Glo)来确定细胞活力。

d2<-fread("GDSC2_fitted_dose_response_25Feb20.csv")%>%as.data.frame()
#GDSC2-LUAD=8897

dd<-subset(d2,d2$TCGA_DESC=="LUAD")

dd$IC50<-exp(log(dd$LN_IC50))
#把ln(IC50)还原成IC50

dd2<-dd[,c(6,9,10,11,16,20)]

###读入细胞系信息

cc<-fread("model_list_20210324.csv")%>%as.data.frame()

cc2<-cc[,c(1,2,25,28,42)]
cc2<-subset(cc2,cc2$model_id%in%dd2$SANGER_MODEL_ID)

cc3<-cc2[!is.na(cc2$age_at_sampling),]
#保留有年龄信息的LUAD细胞系 中位年龄是51

dd3<-subset(dd2,dd2$SANGER_MODEL_ID%in%cc3$model_id)
colnames(dd3)[1]<-"model_id"
#药敏数据与细胞系数据的ID对应

cir_df<-inner_join(cc3,dd3) #药敏数据和细胞系年龄信息整合

cir_df2<-cir_df[!is.na(cir_df$IC50),]

cir_df2$group=ifelse(cir_df2$age_at_sampling>60,"old","young")%>%as.factor()
#按年龄60分两组

drug<-cir_df2$DRUG_NAME%>%unique() ###173种药物数据

#####计算每个药物在两组细胞系中的IC50差异
p_60<-c()
for(i in 1:length(drug)){
  
  f<-subset(cir_df2,cir_df2$DRUG_NAME==drug[i])
  t<-table(f$group)
  
  if(t[2]!=0&t[1]!=0){
    w<-wilcox.test(IC50 ~ group, data = f)
    p_60[i]=w$p.value
    names(p_60)[i]<-drug[i]
  }
}

p60_sig<-p_60[p_60<0.05]
p60_sig<-p60_sig[!is.na(p60_sig)] #20种显著差异药物

######显著差异的药物画IC50-boxplot

setwd("E:/Rwork/lung/R5/drug/boxplot/")

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

plot_list<-list()
for(i in 1:length(p60_sig)){
  
  df<-subset(cir_df2,cir_df2$DRUG_NAME==names(p60_sig)[i])
  plot_list[[i]]<-ggboxplot(df, x = "group", y = "IC50",color = "group",add.params = list(color = "group",size=1),add = "jitter",
                 palette = c("#1F78B4","#A65628"))+stat_compare_means(method = "t.test")+ylab(paste("IC50 of ",names(p60_sig)[i]))+xlab("")+theme(legend.position = "none")
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]])

#######与本文衰老相关的靶基因对应的药物-靶点-通路-plot

res_df<-subset(cir_df2,cir_df2$DRUG_NAME %in% names(p60_sig))
write.csv(res_df,"GDSC2_drug_age60.csv")

gg<-read.table("E:/Rwork/lung/immunosense_gene.txt")
##90个免疫衰老基因
library(stringr)

d1_list<-str_split(d1$PUTATIVE_TARGET,",")%>%unlist()

intersect(d1_list,gg$V1) #GDSC1靶基因与免疫衰老基因交集
#"IGF1R" "MTOR"  "ATM"   "BCL2"  "JAK3"  "JAK2"  "SIRT1" "FAS"  

d2_list<-str_split(d2$PUTATIVE_TARGET,",")%>%unlist()

intersect(d2_list,gg$V1) #GDSC2靶基因与免疫衰老基因交集
#"BCL2"  "ATM"   "IGF1R" "JAK2"  "TP53"  "TERT" 

intersect(str_split(res_df$PUTATIVE_TARGET,",")%>%unlist(),gg$V1)
#IGF1R  显著药物靶点与免疫衰老基因交集

IGF1R<-subset(res_df,res_df$PUTATIVE_TARGET=="IGF1R, IR")
#GSK1904529A-"IGF1R, IR"-IGF1R signaling
write.csv(IGF1R,"IGF1R-GSK1904529A.csv")
