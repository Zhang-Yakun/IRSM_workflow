
# GSE11969 ----------------------------------------------------------------


setwd("D:/workplace/geneset/geo/GSE11969/")

library(data.table)
library(dplyr)


# 生存数据预处理 -----------------------------------------------------------------


sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
age=sur[11,]
status=sur[17,]
time=sur[18,]

sur.df<-rbind(sample,age,status,time)
sur.df<-sur.df[,1:90]
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","age","status","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$age<-sur.df2$age %>%str_remove("Age: ")%>%as.numeric()
sur.df2$status<-sur.df2$status %>%str_remove("Status: ")
sur.df2$time<-sur.df2$time%>%str_remove("Survival \\(days\\)\\: ")%>%as.numeric()

sur.df2$status<-ifelse(sur.df2$status=="Alive",0,ifelse(sur.df2$status=="Dead",1,NA))%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

med<-median(sur.df3$age)

sur.df3$group<-ifelse(sur.df3$age>med,"old","young")
table(sur.df3$group)

# 以年龄分组-生存分析 --------------------------------------------------------------
library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

fit<- survfit(Surv(time, status) ~ group, data = sur.df3)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#00BFC4","#F8766D"))


# 读入表达谱 -------------------------------------------------------------------

exp<-fread("GSE11969_series_matrix.txt")%>%as.data.frame()
#探针表达谱

#####4个风险基因
gene4<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = T)%>%as.data.frame()
gene<-gene4$geneids

#####GPL

gpl<-fread('GPL7015.txt',header=T)%>%as.data.frame()
gpl2<-gpl[,c(1,6)]

gpl3<-subset(gpl2,gpl2$`Gene symbol`%in%gene)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$`Gene symbol`),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$sample,]
identical(sur.df3$sample,rownames(exp6))

# 计算风险得分 ------------------------------------------------------------------

ex<-t(exp6)

final_gene<-rownames(ex)
rownames(gene4)<-gene4$geneids
coef<-gene4[final_gene,]
coef<-coef$index.min

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE11969.csv")

fit<- survfit(Surv(time, OS) ~ risk.group, data = rs.df)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = "aaas")


# ROC一年生存率 ----------------------------------------------------------------

library("pROC")
library("survivalROC")
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

rdat<-data.frame(time=rs.df$time,Group=rs.df$risk.group,data=rs.df)
rdat$Group<-ifelse(rdat$Group=="High_risk",1,0)

nobs<-nrow(rdat)
cutoff<-365

myroc<-survivalROC(Stime = rdat$time,status = rdat$Group, marker = rdat$data.risk.score,
                   predict.time = cutoff,span = 0.25*nobs^(-0.20))
require(ggsci)
library("scales")
pal_nejm("default")(8)
show_col(pal_nejm("default")(8))

plot(myroc$FP,myroc$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=("False positive rate"), ##连接
     ylab="True positive rate",
     main="Time dependent ROC")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

legend(0.6,0.2,c(paste("AUC of risk score=",round(myroc$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 

# GSE26939 ----------------------------------------------------------------


setwd("D:/workplace/geneset/geo/GSE26939/")

library(data.table)
library(dplyr)


# 生存数据预处理 -----------------------------------------------------------------


sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
age=sur[10,]

sample<-c()
status<-c()
time<-c()

s<-c("survival_status(0='alive',1='dead'): 1","survival_status(0='alive',1='dead'): 0")

for(i in 1:dim(sur)[2]){
  
  x=sur[,i]
  sample[i]<-x[1]
  
  for(j in 1:length(x)){
    if(x[j]%in%s){
      status[i]<-x[j]
      time[i]<-x[j-1]
    }
  }
}

sur.df<-rbind(sample,age,status,time)
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","age","status","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$age<-sur.df2$age %>%str_remove("age \\(90\\=greater than or equal to 90\\)\\: ")%>%as.numeric()
sur.df2$status<-sur.df2$status %>%str_remove("survival_status\\(0\\=\\'alive\\'\\,1\\=\\'dead\\'\\)\\: ")%>%as.numeric()
sur.df2$time<-sur.df2$time%>%str_remove("survival\\_months\\: ")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

med<-median(sur.df3$age)

sur.df3$group<-ifelse(sur.df3$age>med,"old","young")
table(sur.df3$group)

# 以年龄分组-生存分析 --------------------------------------------------------------
library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

fit<- survfit(Surv(time, status) ~ group, data = sur.df3)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#00BFC4","#F8766D"))


# 读入表达谱 -------------------------------------------------------------------

exp<-fread("GSE26939_series_matrix.txt")%>%as.data.frame()
#探针表达谱

#####4个风险基因
gene4<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = T)%>%as.data.frame()
gene<-gene4$geneids

#####GPL

gpl<-fread('GPL9053.txt',header=T)%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

gpl3<-subset(gpl2,gpl2$ORF%in%gene)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$ORF),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$sample,]
identical(sur.df3$sample,rownames(exp6))

# 计算风险得分 ------------------------------------------------------------------

final_gene<-gene4$geneids
coef<-gene4$index.min

ex<-t(exp6)

write.csv(ex,"exp_GSE26939.csv")

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE26939.csv")

fit<- survfit(Surv(time, OS) ~ risk.group, data = rs.df)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = "aaas")


# ROC一年生存率 ----------------------------------------------------------------

library("pROC")
library("survivalROC")
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2")

rdat<-data.frame(time=rs.df$time,Group=rs.df$risk.group,data=rs.df)
rdat$Group<-ifelse(rdat$Group=="High_risk",1,0)

nobs<-nrow(rdat)
cutoff<-12

myroc<-survivalROC(Stime = rdat$time,status = rdat$Group, marker = rdat$data.risk.score,
                   predict.time = cutoff,span = 0.25*nobs^(-0.20))
require(ggsci)
library("scales")
pal_nejm("default")(8)
show_col(pal_nejm("default")(8))

plot(myroc$FP,myroc$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=("False positive rate"), ##连接
     ylab="True positive rate",
     main="Time dependent ROC")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

legend(0.6,0.2,c(paste("AUC of risk score=",round(myroc$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)# 

