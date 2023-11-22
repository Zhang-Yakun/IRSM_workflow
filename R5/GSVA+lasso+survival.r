
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

# 读入生存数据 ------------------------------------------------------------------

surdata<-fread("Survival_SupplementalTable_S1_20171025_xena_sp")%>%as.data.frame()

surdata2<-surdata[,c("sample","OS","OS.time")]

cdata<-inner_join(gsva_res2,surdata2)
rownames(cdata)<-cdata$sample
cdata<-cdata[,-2] #把sample列变成行名
cdata<-cdata%>%dplyr::select(OS.time,OS,everything())

cdata$sample<-rownames(cdata)
cdata2<-inner_join(score261,cdata)

# cox分析 -------------------------------------------------------------------

Coxoutput=data.frame()

library(survival)
library(survminer)
library(ggpubr)

cox <- coxph(Surv(OS.time, OS) ~ immunosenescence, data = cdata2)
coxSummary = summary(cox)
Coxoutput=cbind(ID="immunosenescence",HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4])

for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}

#Coxoutput <- arrange(Coxoutput,pvalue)   #按照p值排序 [37] pancancer_cox_result
#Coxout <- filter(Coxoutput,Coxoutput$pvalue<0.05)  #[28] P_value significant

write.csv(Coxoutput,'cox_output1.csv', row.names = F)


# 生存分析 --------------------------------------------------------------------

####### age group survival
library(cowplot)

fit1 <- survfit(Surv(OS.time, OS) ~ group, data = cdata2)

ggsurvplot(fit1,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette =  c("#00BFC4","#F8766D")) #颜色风格


####### gsva-score group survival

m<-median(cdata2$immunosenescence)

cdata2$score_group<-ifelse(cdata2$immunosenescence>m,"high-score","low-score")

fit2 <- survfit(Surv(OS.time, OS) ~ score_group, data = cdata2)

ggsurvplot(fit2,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#E7B800","#2E9FDF"))
 #颜色风格,可选调色板有 "grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons"和"rickandmorty".


# lasso regression --------------------------------------------------------

library(glmnet)
library(survival)

options(stringsAsFactors = FALSE) #禁止chr转成factor

rownames(exp261)<-exp261$V1
ex<-exp261[,intersect(geneset$gene,colnames(exp261))]
ex<-t(ex)
ex<-ex[c("CD40LG","TLR2","IL7","TLR3","SELL","TLR7","CX3CR1"),]
su<-data.frame(time=cdata2$OS.time,OS=cdata2$OS)
rownames(su)<-cdata2$sample  #261
su<-subset(su,su$time!=0) #去掉生存时间为零的样本
ex<-ex[,rownames(su)]

identical(colnames(ex),rownames(su)) ##检查样品信息跟表达量矩阵的样品是否一致

cvfit = cv.glmnet(t(ex), Surv(su$time,su$OS), 
                  #10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提
                  nfold=10,
                  family = "cox") 
plot(cvfit)

fit <- glmnet(t(ex), Surv(su$time,su$OS), 
              family = "cox") 
plot(fit, label = TRUE,xlab="log Lambda")

cvfit$lambda.1se
coef.min = coef(cvfit, s = "lambda.min") 
active.min = which(coef.min != 0)
geneids <- rownames(ex)[active.min]
geneids  #提取选中的基因名
index.min = coef.min[active.min]
index.min  #提取基因对应的系数coef
combine <- cbind(geneids, index.min)
write.csv(combine,"lasso_gene_index.csv",row.names = F)


# 计算风险得分 ------------------------------------------------------------------

final_gene<-geneids
coef<-index.min
expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=su$time,OS=su$OS,risk.score=risk.score,group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_result.csv")

fit<- survfit(Surv(time, OS) ~ group, data = rs.df)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = "aaas")

diffr<-survdiff(Surv(time, OS) ~ group, data = rs.df)

write.csv(rs.df,"risk_score.csv")


# ROC 一年生存率 --------------------------------------------------------------

rdat<-data.frame(time=rs.df$time,Group=rs.df$group,data=rs.df)
rdat$Group<-ifelse(rdat$Group=="High_risk",1,0)

gsva.df<-inner_join(score261,surdata2)
rdat2<-data.frame(time=gsva.df$OS.time,Group=gsva.df$immunosenescence,data=gsva.df)
rdat2$Group<-ifelse(rdat2$Group>median(rdat2$data.immunosenescence),1,0)

library("pROC")
library("survivalROC")
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

nobs<-nrow(rdat)
cutoff<-365

myroc<-survivalROC(Stime = rdat$time,status = rdat$Group, marker = rdat$data.risk.score,
                   predict.time = cutoff,span = 0.25*nobs^(-0.20))

nobs2<-nrow(rdat2)
cutoff<-365
myroc2<-survivalROC(Stime = rdat2$time,status = rdat2$Group, marker = rdat2$data.immunosenescence,
                   predict.time = cutoff,span = 0.25*nobs2^(-0.20))


require(ggsci)
library("scales")
pal_nejm("default")(8)
show_col(pal_nejm("default")(8))

## risk-score
plot(myroc$FP,myroc$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=("False positive rate"), ##连接
     ylab="True positive rate",
     main="Time dependent ROC")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

## gsva-score
lines(myroc2$FP,myroc2$TP, type="l",col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
legend(0.6,0.2,c(paste("AUC of Risk score =",round(myroc$AUC,3)),
                 paste("AUC of GSVA score =",round(myroc2$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 

####

gsva.df2<-data.frame(gsva.score=gsva.df$immunosenescence,sample=gsva.df$sample)
rs.df$sample<-rownames(rs.df)
score2df<-inner_join(rs.df,gsva.df2)

rownames(score2df)<-score2df$sample

surdata_luad<-subset(surdata,surdata$`cancer type abbreviation`=="LUAD")
cli<-dplyr::select(surdata_luad,sample,age_at_initial_pathologic_diagnosis,
                   gender,ajcc_pathologic_tumor_stage,OS,OS.time)

colnames(cli)[2]<-"age"
colnames(cli)[4]<-"tumor_stage" 
rownames(cli)<-cli$sample

cli259<-cli[score2df$sample,]

identical(rownames(cli259),rownames(score2df))
cli259$gsva.score<-score2df$gsva.score
cli259$risk.score<-score2df$risk.score

write.csv(cli259,"gsva_risk_score.csv",row.names = F)

# PCA ---------------------------------------------------------------------

######GSVA-score-PCA

for(i in 1:dim(ex)[1]){
  
  ex[i,is.na(ex[i,])]=0
  
}#去掉表达谱里的NA

group_list<-data.frame(group=gsva.df$group)
#20个疾病样本的聚类分组信息
rownames(group_list)<-gsva.df$sample

pex<-t(ex)
sample259<-intersect(rownames(pex),rownames(group_list))

group_list<-group_list[sample259,]
group_list<-data.frame(group=group_list)
rownames(group_list)<-sample259

pex<-pex[rownames(group_list),]

identical(rownames(group_list),rownames(pex))
#分组信息按表达谱的样本顺序排序

input_expr<-pex
input_meta<-data.frame(subtype=group_list$group)
rownames(input_meta)<-rownames(input_expr)

pca.results <- prcomp(input_expr, center = TRUE, scale. = FALSE)

library(devtools)
library(ggord)
library(ggplot2)
library(plyr)
source('F:/shan_work/cluster/geom_ord_ellipse.R') 

ggord(pca.results,
      grp_in = input_meta$subtype, repel = TRUE,
      ellipse = FALSE,
      size = 2,
      alpha = 0.5,
      cols = c("#223D6C", "#D20A13"),
      arrow = NULL, txt = NULL
) +
  theme(panel.grid = element_blank()) +
  geom_ord_ellipse(
    ellipse_pro = 0.95,
    size = 1,
    lty = 1
  )

######risk-score-PCA

pex2<-pex[,geneids]

identical(rownames(group_list),rownames(pex2))
#分组信息按表达谱的样本顺序排序

input_expr<-pex2
input_meta<-data.frame(subtype=group_list$group)
rownames(input_meta)<-rownames(input_expr)

pca.results <- prcomp(input_expr, center = TRUE, scale. = FALSE)

library(devtools)
library(ggord)
library(ggplot2)
library(plyr)
source('F:/shan_work/cluster/geom_ord_ellipse.R') 

ggord(pca.results,
      grp_in = input_meta$subtype, repel = TRUE,
      ellipse = FALSE,
      size = 2,
      alpha = 0.5,
      cols = c("#223D6C", "#D20A13"),
      arrow = NULL, txt = NULL
) +
  theme(panel.grid = element_blank()) +
  geom_ord_ellipse(
    ellipse_pro = 0.95,
    size = 1,
    lty = 1
  )


# 风险基因热图--------------------------------------------------

rs.df$sample<-rownames(rs.df)
ex2<-t(ex)%>%as.data.frame()
ex2$sample<-rownames(ex2)

risk_df<-inner_join(rs.df,ex2)
risk_df<-arrange(risk_df,group,risk.score)

###
anno_col = data.frame(
  group = as.factor(risk_df$group)
)

rownames(anno_col) = risk_df$sample

exp_risk<-risk_df[,c("TLR2","IL7","CX3CR1","TLR3","CD40LG")]
rownames(exp_risk)<-risk_df$sample

heat_in<-t(exp_risk)

library(pheatmap)
pheatmap(heat_in,cellwidth = 1, cellheight = 18, fontsize = 8,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("navy", "white", "firebrick2"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = anno_col,
         treeheight_col = "0",#不画树
         border_color = "NA")

####

library(ggplot2)

vio_input<-risk_df
vio_input$group<-as.factor(vio_input$group)

t.test(vio_input$CD40LG~vio_input$response)[[3]]
t.test(vio_input$CX3CR1~vio_input$response)[[3]]
t.test(vio_input$IL7~vio_input$response)[[3]]
t.test(vio_input$TLR3~vio_input$response)[[3]]
t.test(vio_input$TLR2~vio_input$response)[[3]]