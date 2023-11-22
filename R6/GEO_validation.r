
setwd("E:/Rwork/lung/GSE72094_LUAD/")
#	GPL15048

library(data.table)
library(dplyr)


# 生存数据预处理 -----------------------------------------------------------------


sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
gender=sur[13,]
age=sur[14,]
status=sur[18,]
time=sur[19,]
                   
sur.df<-rbind(sample,gender,age,status,time)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","gender","age","status","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$gender<-sur.df2$gender%>%str_remove("gender: ")
sur.df2$age<-sur.df2$age %>%str_remove("age_at_diagnosis: ")
sur.df2$status<-sur.df2$status %>%str_remove("vital_status: ")
sur.df2$time<-sur.df2$time %>%str_remove("survival_time_in_days: ")

sur.df2$status<-ifelse(sur.df2$status=="Alive",0,ifelse(sur.df2$status=="Dead",1,NA))
sur.df2$age<-as.numeric(sur.df2$age)

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值
sur.df3$group<-ifelse(sur.df3$age>70,"old","young")
table(sur.df3$group)
sur.df3$time<-as.numeric(sur.df3$time)

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

exp<-fread("GSE72094_series_matrix.txt")%>%as.data.frame()
#探针表达谱

#####12个风险基因
gene12<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = F)
gene<-gene12$V1

#####GPL15048
library(Biobase)
gpl<-fread('GPL15048.txt',header=T)%>%as.data.frame()
gpl2<-gpl[,c(1,4)]

gpl3<-subset(gpl2,gpl2$GeneSymbol%in%gene)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$GeneSymbol),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$sample,]
identical(sur.df3$sample,rownames(exp6))

# 计算风险得分 ------------------------------------------------------------------

coef.df<-fread("lasso_gene_index.csv")%>%as.data.frame()
coef.df<-subset(coef.df,coef.df$geneids%in%colnames(exp6))

final_gene<-coef.df$geneids
coef<-coef.df$index.min

ex<-t(exp6)

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE72094.csv")

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


# GEO68465 ----------------------------------------------------------------

setwd("E:/Rwork/lung/GSE68465_LUAD/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
gender=sur[10,]
age=sur[11,]
status=sur[13,]
time=sur[20,]

sur.df<-rbind(sample,gender,age,status,time)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","gender","age","status","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$gender<-sur.df2$gender%>%str_remove("Sex: ")
sur.df2$age<-sur.df2$age %>%str_remove("age: ")%>%as.numeric()
sur.df2$status<-sur.df2$status %>%str_remove("vital_status: ")
sur.df2$time<-sur.df2$time %>%str_remove("months_to_last_contact_or_death: ")%>%as.numeric()

sur.df2$status<-ifelse(sur.df2$status=="Alive",0,ifelse(sur.df2$status=="Dead",1,NA))
sur.df2$age<-as.numeric(sur.df2$age)

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

median(sur.df3$age)
sur.df3$group<-ifelse(sur.df3$age>70,"old","young")
table(sur.df3$group)
sur.df3$time<-as.numeric(sur.df3$time)

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


######表达谱

exp<-fread("GSE68465_series_matrix.txt")%>%as.data.frame()
#探针表达谱

gpl<-fread("GPL96.annot")%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

gpl3<-subset(gpl2,gpl2$`Gene symbol` %in%gene)

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

coef.df<-fread("E:/Rwork/lung/GSE72094_LUAD/lasso_gene_index.csv")%>%as.data.frame()
coef.df<-subset(coef.df,coef.df$geneids%in%colnames(exp6))

final_gene<-coef.df$geneids
coef<-coef.df$index.min

ex<-t(exp6)

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE68465.csv")

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

# GSE68571 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE68571/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=colnames(sur)
gender=sur[9,]
age=sur[10,]
status=sur[19,]
time=sur[18,]

sur.df<-rbind(sample,gender,age,status,time)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","gender","age","status","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$gender<-sur.df2$gender%>%str_remove("Sex: ")

sur.df2<-subset(sur.df2,sur.df2$gender!="--")

sur.df2$age<-str_remove(sur.df2$age,pattern="[age (years):]*")%>%as.numeric()
sur.df2$status<-str_sub(sur.df2$status,25,25)%>%as.numeric()
sur.df2$time<-sur.df2$time %>%str_remove("[followup time (months):]*")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

med<-median(sur.df3$age)
sur.df3$group<-ifelse(sur.df3$age>med,"old","young")
table(sur.df3$group)
sur.df3$time<-as.numeric(sur.df3$time)

sur.df3<-sur.df3[,-1]
sur.df3$group<-as.factor(sur.df3$group)

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


######表达谱


#####4个风险基因
gene12<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = T)%>%as.data.frame()
gene<-gene12$geneids
#####
exp<-fread("GSE68571_series_matrix.txt")%>%as.data.frame()
#探针表达谱

gpl<-fread("GPL80.annot")%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

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
exp6<-exp6[rownames(sur.df3),]
identical(rownames(sur.df3),rownames(exp6))

# 计算风险得分 ------------------------------------------------------------------

coef.df<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = T)%>%as.data.frame()

coef.df<-subset(coef.df,coef.df$geneids%in%colnames(exp6))

final_gene<-coef.df$geneids
coef<-coef.df$index.min

ex<-t(exp6)

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE68571.csv")

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

# GSE31852 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE31852/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample<-c()
status<-c()
time<-c()

s<-c("pfsc (1=progressed; 0=not progressed): 1","pfsc (1=progressed; 0=not progressed): 0","progression-free survival status: 1","progression-free survival status: 0")

for(i in 1:dim(sur)[2]){
  
  x=sur[,i]
  sample[i]<-x[1]
  
  for(j in 1:length(x)){
    if(x[j]%in%s){
     status[i]<-x[j]
     time[i]<-x[j+1]
    }
  }
}


sur.df<-rbind(sample,status,time)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","status","time")
sur.df2<-as.data.frame(sur.df2)

sur.df2<-subset(sur.df2,!is.na(sur.df2$status))

status2<-str_split(sur.df2$status,pattern = ":")%>%unlist()
sur.df2$status<-status2[seq(0,length(status2),2)]

time2<-str_split(sur.df2$time,pattern = ":")%>%unlist()
sur.df2$time<-time2[seq(0,length(time2),2)]%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值


######表达谱

#####
exp<-fread("GSE31852_series_matrix.txt")%>%as.data.frame()
#探针表达谱

gpl<-fread("GPL6244.annot")%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

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
exp6<-exp6[rownames(sur.df3),]
identical(rownames(sur.df3),rownames(exp6))

# 计算风险得分 ------------------------------------------------------------------

coef.df<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = T)%>%as.data.frame()
coef.df<-subset(coef.df,coef.df$geneids%in%colnames(exp6))

final_gene<-coef.df$geneids
coef<-coef.df$index.min

ex<-t(exp6)

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE31852.csv")

rs.df$OS<-as.numeric(rs.df$OS)

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

# GEO87340 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE87340/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
gender=sur[11,]
age=sur[12,]
status=sur[17,]
time=sur[18,]

sur.df<-rbind(sample,gender,age,status,time)
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","gender","age","status","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$gender<-sur.df2$gender%>%str_remove("Sex: ")
sur.df2$age<-sur.df2$age %>%str_remove("age: ")%>%as.numeric()
sur.df2$status<-sur.df2$status %>%str_remove("survstatus: ")%>%as.numeric()
sur.df2$time<-sur.df2$time %>%str_remove("time to last followup: ")%>%as.numeric()
###

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

med<-median(sur.df3$age)
sur.df3$group<-ifelse(sur.df3$age>med,"old","young")
table(sur.df3$group)
sur.df3$time<-as.numeric(sur.df3$time)

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


######表达谱
gene12<-fread("D:/workplace/geneset/geo/lasso_gene_index.csv",header = T)%>%as.data.frame()
gene<-gene12$geneids
###
exp<-fread("GSE87340_RPKM_log2.txt")%>%as.data.frame()
#探针表达谱

exp2<-subset(exp,exp$ID%in%gene)

rownames(sur.df3)<-str_remove(rownames(sur.df3),pattern = ".RNA-seq")

rownames(exp2)<-exp2$ID
exp2<-exp2[,-1]
exp6<-t(exp2)%>%as.data.frame()
exp6<-exp6[rownames(sur.df3),]
identical(rownames(sur.df3),rownames(exp6))

# 计算风险得分 ------------------------------------------------------------------

final_gene<-gene12$geneids
coef<-gene12$index.min

ex<-t(exp6)

expr_train<-ex[final_gene,]
risk.score<-c()
for(i in 1:ncol(expr_train)){
  risk.score[i]<-sum(coef*expr_train[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(time=sur.df3$time,OS=sur.df3$status,risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))
rownames(rs.df)<-colnames(expr_train)
write.csv(rs.df,"RiskScore_GSE87340.csv")

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
