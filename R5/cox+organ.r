

# 单因素cox-临床信息 -------------------------------------------------------------

cli259<-read.csv("gsva_risk_score.csv")
rownames(cli259)<-cli259$sample
cli259<-cli259[,-1]

cli9<-cli259 %>% dplyr::select(OS.time,OS,everything())

cli9$tumor_stage<-ifelse(cli9$tumor_stage %in% c("Stage IIIA","Stage IIIB","Stage IV"),1,0)
#stage3,4=1（严重），stage1,2=0(轻微)

cop=data.frame()

for(i in colnames(cli9[,3:ncol(cli9)])){
  cox <- coxph(Surv(OS.time, OS) ~ cli9[,i], data = cli9)
  coxSummary = summary(cox)
  cop=rbind(cop,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                      z=coxSummary$coefficients[,"z"],
                      pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                      lower=coxSummary$conf.int[,3],
                      upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  cop[,i] <- as.numeric(as.vector(cop[,i]))
}

cop <- arrange(cop,desc(pvalue))   #按照p值排序

write.csv(cop,'clinical_cox_output.csv', row.names = F)

tabletext<-cbind(c("",cop$gene),
                 c("pvalue",round(cop$pvalue,3)),
                 c("Hazard Ratio",paste(format(round(cop$HR,3),nsmall=2)," (",format(round(cop$lower,3),nsmall = 2)," - ",format(round(cop$upper,3),nsmall = 2),")",sep="")))
head(tabletext)
library(forestplot)
forestplot(labeltext=tabletext,mean=c(NA,round(cop$HR,3)),
           lower=c(NA,round(cop$lower,3)),upper=c(NA,round(cop$upper,3)),col=fpColors(box="green", lines="steelblue", zero = "black"),#box颜色
           xlab= "Hazard Ratio",zero=1,boxsize=0.2,lwd.ci=2 ,#boxsize=c(NA,NA,Coxoutput$HR,NA)/10
)

# 临床信息-多因素cox回归分析 ----------------------------------------------------------

library(survival)
library(survminer)

s=Surv(OS.time,OS)~gender+tumor_stage+gsva.score+risk.score
model<-survival::coxph(s,data = cli9)
options(scipen = 1)
ggforest(model, data = cli9, main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7,
         refLabel = "reference", noDigits = 2)

model<-coxph(Surv(OS.time,OS)~gender+tumor_stage+gsva.score+risk.score,data=cli9)
ggforest(model)


# organ 组织分布图 -------------------------------------------------------------

library(gridExtra)
library(ggpolypath)
library(ggplot2)

devtools::install_github("jespermaag/gganatogram")
library(gganatogram)

or<-fread("CSNK1G3_tissues_MERAV.txt")%>%as.data.frame()
or$organ<-str_extract(or$GENE,"\\b[a-zA-Z]+")

s<-str_split(or$GENE,"_")
exam<-c()
exam2<-c()
for(i in 1:nrow(or)){
  exam[i]=s[[i]][2]
  exam2[i]=s[[i]][1]
}
or$organ<-exam2
or$type<-exam
or$organ= str_replace_all(or$organ,"[.]","_")
or$organ= str_to_lower(or$organ)

myorgan<-intersect( hgMale_key$organ,or$organ)%>%unique()
myor<-filter(or,or$organ%in%myorgan)
myor<-arrange(myor,myor$organ)
myt<-table(myor$organ) #14 organ
mycol<-rep(rainbow(14),as.numeric(myt))
myor$color<-mycol
organ<-myor[,2:5]
organ$type<-str_sub(organ$type,1,13)
organ$value<-organ$CSNK1G3

hgmale <- gganatogram(data=organ, 
                      fillOutline='white', #没有organ的位置用白色填充
                      organism='human', 
                      sex='male', #性别
                      fill='value') + #还可以用color列填充
  facet_wrap(~as.factor(type)) + #对比tumor和normal
  scale_fill_gradient(low = "black", high = "red") + 
  labs(fill = "value") + 
  theme_void() #不画坐标轴

write.csv(organ[,2:5],"CSNK1G3_organ.csv",row.names = F)
