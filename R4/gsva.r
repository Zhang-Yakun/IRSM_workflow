
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)

setwd("E:/Rwork/lung/R5/")

# mygeneset ---------------------------------------------------------------

geneset<-clusterProfiler::read.gmt("E:/Rwork/lung/my_geneset.gmt")

mypathway<-list()

mypathway[[1]]<-geneset$gene
names(mypathway)[1]<-"immunosenescence"


# 读入基因表达数据 ----------------------------------------------------------------

tpm<-fread("LUAD_tpm_age.csv")%>%as.data.frame()
#所有基因表达log2+1(TPM值)

group<-data.frame(sample=tpm$sample,age=tpm$age)

tpm_exp<-tpm[,3:17559]
rownames(tpm_exp)<-tpm$sample
tpm_exp<-t(tpm_exp)

# GSVA分析 ------------------------------------------------------------------
gsva_res<- gsva(as.matrix(tpm_exp), mypathway)
head(gsva_res)  #GSVA result

gsva_res2<-t(gsva_res)%>%as.data.frame()
gsva_res2$sample <- rownames(gsva_res2)


# 差异分析-GSVA score在组间的差异 ------------------------------------------------

###样本分组
exp261<-fread("exp_261.csv")%>%as.data.frame()

sample_group<-data.frame(sample=exp261$V1,group=exp261$group)
score261<-inner_join(sample_group,gsva_res2)

wilcox.test(immunosenescence ~ group, data = score261)
#p-value = 0.0002871

write.csv(score261,"GSVA_score.csv",row.names = F)

###histgram-GSVA-score

library(ggplot2)
library(reshape2)

#输入包含2个维度的数据。每一行是一簇柱子，每一列用不同的颜色图例区分

rownames(score261)<-score261$sample
df_gsva_score<-score261[,-1]

df_gsva_score$group<-as.factor(df_gsva_score$group)

df_gsva_score<-arrange(df_gsva_score,group,immunosenescence)

box_young<-subset(df_gsva_score,df_gsva_score$group=="young")
box_old<-subset(df_gsva_score,df_gsva_score$group=="old")


##背靠背直方图

h1=hist(box_old$immunosenescence,plot = F)
h2=hist(box_young$immunosenescence,plot = F)
#绘制两个直方图,数据存在h1h2两个对象中
h2$counts= - h2$counts
#将h2的值反过来
hmax= max(h1$counts)
hmin= min(h2$counts)
#设置y轴取值范围
X = c(h1$breaks,h2$breaks)
xmax = max(X)
xmin = min(X)
#设置x轴取值范围
plot(h1,ylim = c(hmin,hmax),col="#00AFBB",
     xlim = c(xmin,xmax),main = 'GSVA score of immunosenescence gene set')
lines(h2,col='#FC4E07')

#####年龄横轴，GSVA-score纵轴的相关性散点图

age_df<-data.frame(sample=exp261$V1,age=exp261$age)

gsva_age<-inner_join(age_df,score261)
rownames(gsva_age)<-gsva_age$sample

library(ggplot2)

theme <- theme(panel.background = element_blank(), # 去掉背景格子
               # 显示x平行网格线
               #panel.grid.major.x = element_line(colour = "black"), 
               # 显示x轴坐标
               axis.line.x = element_line(colour = "black"),
               axis.line.y = element_line(colour = "black"),
               axis.title.y = element_blank())

b <- ggplot(gsva_age, aes(x = age, y = immunosenescence))
# Scatter plot with regression line

b + geom_point(aes(color=group))+
  geom_smooth(aes(color = group, fill = group), method = "lm",fullrange = TRUE)+
  facet_wrap(~group) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme_bw()

cor<-cor.test(as.numeric(gsva_age$age), as.numeric(gsva_age$immunosenescence) , method="spearman")
cor$p.value
cor$estimate

######gsva_score_corplot2
b + geom_point(aes(color = group, shape = group)) +
  geom_rug(aes(color =group)) +
  geom_smooth(aes(color = group), method = lm, 
              se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#00AFBB",  "#FC4E07"))+
  ggpubr::stat_cor(aes(color = group), label.x = 60)+theme
