
setwd("E:/Rwork/lung/R5/")
options(stringsAsFactors = F)

library(oncoPredict)
library(gtools)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
# 用system.time来返回计算所需时间
#system.time({ })

#样本分组和表达谱

exp261<-fread("exp_261.csv")%>%as.data.frame()

rownames(exp261)<-exp261$V1
exp<-exp261[,5:17561]

group<-fread("risk_score.csv")%>%as.data.frame()

exp2<-exp[group$V1,]
Expr<-as.matrix(exp2)

Expr2<-t(Expr)%>%as.matrix()
#以GDSC2-训练集的表达量矩阵和药物处理信息
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

GDSC2_Expr = readRDS(file="GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res = readRDS(file = "GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 

#预测药物反应
  
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = Expr2,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

####分组比较药物差异
drug<-read.csv("E:/Rwork/lung/R5/calcPhenotype_Output/myDrugPredictions.csv")
library(dplyr)
library(data.table)
rownames(drug)<-drug$X
drug<-drug[,-1]

group<-fread("E:/Rwork/lung/R5/risk_score.csv")%>%as.data.frame()
drug$group<-group$group

drug2<-drug%>%dplyr::select(group,everything())

library(stringr)
colnames(drug2)<-str_remove(colnames(drug2),"_[0-9]*")
colnames(drug2)<-str_replace(colnames(drug2),"[._]","-")
###

pt<-c()

for(i in 2:dim(drug2)[2]){
  
  f<-drug2[,i]

  t<-t.test(f ~ drug2$group)
  pt[i-1]=t$p.value
  names(pt)[i-1]<-colnames(drug2)[i]
    
}

p_sig<-p[p<0.01]
p_sig<-p_sig[!is.na(p_sig)]

pt_sig<-pt[pt<0.01]

names(pt_sig)

#p_sig<-p_sig[-30]
#p_sig<-p_sig[-30]
#p_sig<-p_sig[-30]

write.csv(drug2,"E:/Rwork/lung/R5/drug/oncopredict/drug_group.csv")
res<-as.data.frame(pt_sig)
res$drug<-names(pt_sig)
write.table(res,"E:/Rwork/lung/R5/drug/oncopredict/drug_pvalue38.txt",row.names = F,quote = F)

######显著差异的药物画IC50-boxplot

setwd("E:/Rwork/lung/R5/boxplot/")

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

plot_list<-list()
for(i in 1:length(pt_sig)){
  
  df<-drug2[,names(pt_sig)[i]]%>%as.data.frame()
  colnames(df)<-names(pt_sig)[i]
  
  df2<-cbind(drug2[,1],df)
  colnames(df2)[1]<-"group"
  plot_list[[i]]<-ggboxplot(df2, x = "group", y = names(pt_sig)[i] , color = "group",add.params = list(color = "group",size=1),add = "jitter",
                            palette = c("#1F78B4","#A65628"),title=names(pt_sig)[i])+stat_compare_means(method = "t.test")+ylab("Estimated IC50")+xlab("")+theme(legend.position = "none")
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]])

plot_grid(plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],plot_list[[26]],
          plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],
          plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]],plot_list[[37]],plot_list[[38]])

# sankeyplot --------------------------------------------------------------

gdsc<-fread("E:/Rwork/lung/R5/drug/GDSC2_fitted_dose_response_25Feb20.csv")%>%as.data.frame()

luad<-subset(gdsc,gdsc$TCGA_DESC=="LUAD")

drug_name<-intersect(luad$DRUG_NAME,names(pt_sig))

san_input<-subset(luad,luad$DRUG_NAME%in%names(pt_sig))

san_df<-data.frame(drug=san_input$DRUG_NAME,
                   target=san_input$PUTATIVE_TARGET,
                   pathway=san_input$PATHWAY_NAME)

library(ggplot2)
library(RColorBrewer)
library(ggalluvial)

display.brewer.all()
brewer.pal.info

col<- c(c(brewer.pal(9,"YlOrRd"),"#662506"),'#2e1f54', '#52057f',
        "#99D8C9","#66C2A4","#41AE76","#238B45","#006D2C",
        "#00441B","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5",
        "#08519C","#08306B","#3288BD","#5E4FA2","#3288BD","#66C2A5","#66C2A5")#自定义颜色
#新建一个PDF文件
pdf("test.pdf",width = 8, height = 6)
#绘图

#格式转换——通过to_lodes_form函数将数据转换为作图所需要的数据 
df <- to_lodes_form(san_df[,1:ncol(san_df)],
                    axes = 1:ncol(san_df),
                    id = "value")
print(df)#预览数据

ggplot(df, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.5,#透明度
            color = 'white',#间隔颜色
            size = 0.1)+#间隔宽度
  geom_stratum(width = 0.28)+#图中方块的宽度
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  #scale_fill_manual(values = col)+#自定义颜色
  theme_void()+#主题（无轴及网格线）
  theme(legend.position = 'none')#去除图例
dev.off()#关闭PDF 

# LAG3, PDCD1, CTLA-4,TIGIT在两组中的差异boxplot ------------------------------------------------------

e<-read.csv("E:/Rwork/lung/R5/exp_261.csv")

rownames(e)<-e$X
e<-e[,-1]

e2<-e[,c("LAG3","PDCD1","CTLA4","TIGIT")]

group<-fread("E:/Rwork/lung/R5/risk_score.csv")%>%as.data.frame()

e3<-e2[group$V1,]
identical(group$V1,rownames(e3))

e3$group<-group$group

pi<-c()

for(i in 1:4){
  
  f<-e3[,i]
  
  t<-t.test(f ~ e3$group)
  pi[i]=t$p.value
  names(pi)[i]<-colnames(e3)[i]
  
}

#####小提琴图
library(vioplot)
library(data.table)
library(dplyr)
library(stringr)

library(patchwork)

e3$group<-as.factor(e3$group)

par(mfrow=c(1,4))
vioplot(LAG3~group,
        data = e3,
        col=c("#FFD121","#5D90BA"),
        xlab=pi[1], ylab=colnames(e3)[1],
)
vioplot(PDCD1~group,
        data = e3,
        col=c("#FFD121","#5D90BA"),
        xlab=pi[2], ylab=colnames(e3)[2])
vioplot(CTLA4~group,
        data = e3,
        col=c("#FFD121","#5D90BA"),
        xlab=pi[3], ylab=colnames(e3)[3])
vioplot(TIGIT~group,
        data = e3,
        col=c("#FFD121","#5D90BA"),
        xlab=pi[4], ylab=colnames(e3)[4])

#####
library(vioplot)

High=129                                               
Low=130                                                    

library(data.table)
library(dplyr)

rt<-e3[,-5]
#行是样本 列是细胞的免疫细胞丰度结果

h<-subset(e3,e3$group=="High_risk")
l<-subset(e3,e3$group=="Low_risk")

df1<-rt[rownames(h),]
df2<-rt[rownames(l),]

rt2<-rbind(df1,df2)
rt2$group<-c(rep("High",129),rep("Low",130))

df=data.frame()

for(i in 1:4){
  
  x1<-rt2[,i]
  d1<-data.frame(group=rt2$group,gene=rep(colnames(rt2)[i],259),value=x1)
  df<-rbind(df,d1)
}


####

e <- ggplot(df, aes(x = gene, y = value))

e +
  geom_boxplot(
    aes(color = group), width =0.6,
    position = position_dodge(0.9)
  ) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  theme_bw()
