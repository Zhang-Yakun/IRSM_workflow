setwd("E:/Rwork/lung/")

library(data.table)
library(dplyr)

####读取表达和差异结果
exp261<-fread("E:/Rwork/lung/exp_261.csv")%>%as.data.frame()
DEG<-read.csv("E:/Rwork/lung/limma_result_old.csv")

d<-subset(DEG,DEG$adj.P.Val <0.05)

immune_gene90<-read.table("immunosense_gene.txt",header=F) #90个免疫衰老基因列表
immune_gene90<-immune_gene90$V1%>%unique()%>%as.character()

dd<-subset(d,d$X%in%immune_gene90)
dd<-subset(dd,dd$logFC>0)
#老年组中显著上调的免疫衰老基因数：31个(p<0.05), 11个(adj.p<0.05)
gene_11<-dd$X

write.table(gene_11,"老年组上调基因11.txt",row.names = F,col.names = F,quote = F)

deg_exp<-exp261[,intersect(dd$X,colnames(exp261))]
rownames(deg_exp)<-exp261$X

deg_exp$X<-rownames(deg_exp)
group_list<-exp261[,1:2]

deg_exp2<-inner_join(group_list,deg_exp)
deg_exp3<-arrange(deg_exp2,desc(group))
rownames(deg_exp3)<-deg_exp3$X

group_sort<-deg_exp3[,1:2]

deg_exp4<-deg_exp3[,-1]
deg_exp4<-deg_exp4[,-1]

exp_heat<-t(deg_exp4)

identical(group_sort$X,colnames(exp_heat))

####热图展示老年组显著上调的基因
library(pheatmap)
cormat<-round(cor(exp_heat,method = "spearman"),2)

table(group_sort$group)

annotation_col<-data.frame(
  Type = factor(group_sort$group)) #按年龄分组
rownames(annotation_col) = colnames(exp_heat)

ann_colors=list(Type=c(young="#DDA678",old="#88BA9F"))

n <- t(scale(t(exp_heat)))#计算Z-score
n[n > 2] <- 2
n[n < -2] <- -2
df <- n

pheatmap(df,cellwidth = 1, cellheight = 6, fontsize = 6,
         method="pearson", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         #scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("#3178A3","white","#D86472"))(100),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors=ann_colors,
         treeheight_col = "0",#不画树
         border_color = "NA")


# 小提琴图展示老年组和年轻组免疫衰老基因的差异 --------------------------------------------------

vio_exp<-exp261[,intersect(immune_gene90,gene_11)]

vio_exp$mean<-apply(vio_exp,1,mean)
rownames(vio_exp)<-exp261$V1
vio_exp$group<-exp261$group

library(ggstatsplot)
library(ggplot2)
library(ggthemes)
ggbetweenstats(data=vio_exp,x=group,y=mean,title = "Expression of up-regulated immunosenescence genes",
               ggtheme = theme_few(),package="wesanderson",palette="Darjeeling1")


# PCA展示年轻组和老年组的差异 ---------------------------------------------------------

for(i in 1:dim(vio_exp)[1]){
  
  vio_exp[i,is.na(vio_exp[i,])]=0
  
}#去掉表达谱里的NA

group_list<-data.frame(group=vio_exp$group)
#20个疾病样本的聚类分组信息
rownames(group_list)<-rownames(vio_exp)

identical(rownames(group_list),rownames(vio_exp))
#分组信息按表达谱的样本顺序排序

input_expr<-vio_exp[,1:11]
input_meta<-data.frame(subtype=group_list$group)
rownames(input_meta)<-rownames(input_expr)

pca.results <- prcomp(input_expr, center = TRUE, scale. = FALSE)

#install.packages("devtools", type = "win.binary")

library(devtools)

#安装ggord
#install_github('fawda123/ggord')

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

# 免疫浸润分析 ------------------------------------------------------------------

im<-exp261[,-3]
im<-im[,-3]
im_input<-arrange(im,desc(group))
rownames(im_input)<-im_input$V1

im_input2<-t(im_input)%>%as.data.frame()
im_input2<-im_input2[-1,]

im_input2$Symbol<-rownames(im_input2)
im_input3<-im_input2 %>% dplyr::select("Symbol",everything())

write.table(im_input3,"ImmuCellAI_input.txt",row.names = F,quote = F,sep = "\t")
#免疫浸润输入文件

#####violin plot

library(vioplot)
setwd("E:/Rwork/lung/Immucell/")      
young=136                                                
old=125                                                    

library(data.table)
library(dplyr)

rt<-fread("ImmuCellAI_abundance_result.txt")%>%as.data.frame()
#行是样本 列是细胞的免疫细胞丰度结果

rownames(rt)<-rt$V1
rt<-rt[,-1]
rt<-rt[,-25]

#pdf("vioplot.pdf",height=8,width=15)            
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,71),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=20,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，低TMB用绿色表示，高TMB用红色表示
for(i in 1:ncol(rt)){                 #for 1 到 22 列
  if(sd(rt[1:young,i])==0){          #如果从1到96行的1到22列的标准差等于0，rt的第i列第1行等于0.001
    rt[1,i]=0.001                        #
  }
  if(sd(rt[(young+1):(young+old),i])==0){
    rt[(young+1),i]=0.001
  }
  youngData=rt[1:young,i]
  oldData=rt[(young+1):(young+old),i]
  vioplot(youngData,at=3*(i-1),lty=1,add = T,col = '#FFB90F',side = "left",pchMed = 20)
  vioplot(oldData,at=3*(i-1),lty=1,add = T,col = '#008B8B', colMed = 'green',side = "right",pchMed = 18)
  wilcoxTest=wilcox.test(youngData,oldData)
  p=round(wilcoxTest$p.value,2)
  mx=max(c(youngData,oldData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.7)
  text(seq(1,72,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.8,srt = 45,pos=2)
}
#左边黄色-young,右边绿色-old
dev.off()

#############堆叠条形图

library(ggplot2)
library(reshape2)

#输入包含2个维度的数据。每一行是一簇柱子，每一列用不同的颜色图例区分

rt$sample<-rownames(rt)

group_list<-data.frame(sample=rownames(group_list),group=group_list$group)

df_g<-inner_join(group_list,rt)
df_g<-df_g[,-1]

df = reshape2::melt(df_g)
# 把数据转换成ggplot常用的类型（长数据）

p = ggplot(df_g, aes(x=factor(group,levels =unique(group)),  # 转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
                   y=value, 
                   fill=factor(variable,levels = unique(variable)), 
))+
  labs(
    x="",   # 调整x轴名称
    y="",   # 调整y轴名称
    fill="" # 调整图例名称
  )
p + geom_bar(
  position="fill",
  stat="identity"
)+theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ylab("Immune cells infiltration ratio(%)")

# 免疫细胞相关性图 ----------------------------------------------------------------

library(corrplot)
cor_input<-rt[,-25]
#免疫细胞矩阵，列是免疫细胞，行是样本

m<-cor(cor_input)

col3 <- colorRampPalette(c("blue", "white", "red"))#设置颜色
res1 <- cor.mtest(m, conf.level = .95) 

corrplot(corr=m,method = "color",order="AOE",addCoef.col = "black",tl.cex=0.8,number.cex=0.6)

write.csv(m,"免疫细胞比例相关性矩阵.csv")

# DEG与免疫细胞的相关性 ----------------------------------------------------------

options(stringsAsFactors = FALSE) #禁止chr转成factor

rownames(im)<-im$V1
exp_matrix<-im[,-1]
#完整的表达谱数据

gene_11<-read.table("老年组上调基因11.txt")
exp_matrix<-exp_matrix[,gene_11$V1]
exp_matrix<-t(exp_matrix)%>%as.data.frame()

cor_input<-cor_input[colnames(exp_matrix),]
identical(rownames(cor_input),colnames(exp_matrix))

gene <-gene_11$V1
immuscore <- function(gene){
  y <- as.numeric(exp_matrix[gene,])
  colnames <- colnames(cor_input)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(cor_input[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
##先写一个函数,输入一个基因，即可返回跟免疫基因的相关性、p值

immuscore("TLR3")
#以TNF为例，测试一下函数

res <- do.call(rbind,lapply(gene,immuscore))
head(res)
#基因与免疫细胞的相关性结果

#保存到文件
write.csv(res, "cell_gene_correlation.csv", quote = F, row.names = F)

######相关性画图
res$pstar <- ifelse(res$p.value < 0.05,
                    ifelse(res$p.value < 0.01,"**","*"),
                    "")
res$pstar[1:20]

ggplot(res, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size=10,angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 10))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))


# 免疫浸润得分的直方图 ----------------------------------------------------------

box<-data.frame(sample=rownames(rt),score=rt$InfiltrationScore)

box_input<-inner_join(group_sort,box)

box_young<-subset(box_input,box_input$group=="young")
box_young<-arrange(box_young,score)
box_old<-subset(box_input,box_input$group=="old")
box_old<-arrange(box_old,score)

##背靠背直方图

h1=hist(box_old$score,plot = F)
h2=hist(box_young$score,plot = F)
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
plot(h1,ylim = c(hmin,hmax),col="#008B8B",
     xlim = c(xmin,xmax),main = 'Infiltration Score')
lines(h2,col='#FFB90F')


# 三类免疫因子和衰老上调基因的相关性图 ------------------------------------------------------

#####
rownames(im)<-im$V1
exp_matrix<-im[,-1]
exp_matrix<-exp_matrix[,-1]
exp_matrix<-t(exp_matrix)%>%as.data.frame()
#完整的表达谱数据

####免疫因子的表达
Immunoinhibitors<-read.table("Immunostimulators.txt")

vv<-intersect(Immunoinhibitors$V1,rownames(exp_matrix))
#n=23
exp1<-exp_matrix[vv,]
cor_input1<-t(exp1)%>%as.data.frame()

########上调基因的表达

gene_11<-read.table("老年组上调基因11.txt")

gene <-gene_11$V1
immufactor <- function(gene){
  y <- as.numeric(exp_matrix[gene,])
  colnames <- colnames(cor_input1)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(cor_input1[,x]), y , method="spearman")
    data.frame(gene=gene,factor=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
##先写一个函数,输入一个基因，即可返回跟免疫基因的相关性、p值

immufactor("TNF")
#以TNF为例，测试一下函数

res1 <- do.call(rbind,lapply(gene,immufactor))
head(res1)
#基因与免疫因子的相关性结果

#保存到文件
#write.csv(res1, "Chemokines_gene_cor.csv", quote = F, row.names = F)

######相关性画图
res1$pstar <- ifelse(res1$p.value < 0.05,
                     ifelse(res1$p.value < 0.01,"**","*"),
                     "")
res1$pstar[1:20]

ggplot(res1, aes(factor, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size=10,angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 10))+#调整y轴文字
  #调整legend
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"),
       title="Chemokines")


# 挑三组相关性最高的因子图 -----------------------------------------------------------------

c_res1<-subset(res1,res1$p.value<0.05)
c_res1<-subset(c_res1,!c_res1$cor==1)
factor1<-c_res1$factor[which(c_res1$cor==max(abs(c_res1$cor)))]
gene1<-c_res1$gene[which(c_res1$cor==max(abs(c_res1$cor)))]
#最相关的细胞因子和基因


###画相关性图


library(ggstatsplot)#加载包
library(ggside)
p1<-ggscatterstats(data = t(exp_matrix)%>%as.data.frame(),
                   y = TLR7, #gene1
                   x = CD86, #factor1
                   centrality.para = "mean",
                   margins = "both",
                   xfill = "#CC79A7",
                   yfill = "#009E73",
                   marginal.type = "histogram")
p1


# 衰老基因集与特定免疫细胞的相关性散点图 ---------------------------------------------------

####表达数据
exp261<-fread("E:/Rwork/lung/exp_261.csv")%>%as.data.frame()
immune_gene90<-read.table("E:/Rwork/lung/immunosense_gene.txt",header=F) #90个免疫衰老基因列表
immune_gene90<-immune_gene90$V1%>%unique()%>%as.character()

gene_exp<-exp261[,intersect(immune_gene90,colnames(exp261))]
rownames(gene_exp)<-exp261$V1

gene_exp$mean<-apply(gene_exp,1,mean)

####样本分组
group_df<-data.frame(sample=exp261$V1,group=exp261$group)

####免疫细胞表达
rt<-fread("ImmuCellAI_abundance_result.txt")%>%as.data.frame()
#行是样本 列是细胞的免疫细胞丰度结果

rownames(rt)<-rt$V1
rt<-rt[,-1]
rt<-rt[,-25]

rt2<-rt[rownames(gene_exp),]
identical(rownames(gene_exp),rownames(rt2))
geneset_mean<-data.frame(mean=gene_exp$mean)
rownames(geneset_mean)<-rownames(gene_exp)

cor_point<-cbind(geneset_mean,rt2)

cor_df<-data.frame()

for(i in 1:24){
  
  r<-cor.test(cor_point$mean,cor_point[,i+1],method = "pearson")
  cor_r<-r$estimate
  cor_p<-r$p.value

  cor_df2<-data.frame(cor=cor_r,pvalue=cor_p,cell=colnames(cor_point)[i+1])
  cor_df<-rbind(cor_df,cor_df2)
}

cor_dif<-subset(cor_df,cor_df$pvalue<0.05)
cor_dif<-subset(cor_dif,abs(cor_dif$cor)>0.6)
##与免疫衰老基因集表达显著相关的细胞，p<0.05,|r|>0.6
cor_dif$cell

#####相关性散点图
library(ggplot2)
library(ggpubr)
library(ggpmisc)

b1 <- ggplot(cor_point, aes(x = CD4_naive, y = mean))

b1 + geom_point(color="navy")+
  geom_smooth(method = "lm", color = "navy", fill = "lightgray")+
  stat_cor(method = "pearson",label.x = 0, label.y = 4.5) 

b2 <- ggplot(cor_point, aes(x = iTreg, y = mean))

b2 + geom_point(color="navy")+
  geom_smooth(method = "lm", color = "navy", fill = "lightgray")+
  stat_cor(method = "pearson",label.x = 0, label.y = 4.5) 

b3 <- ggplot(cor_point, aes(x = Tfh, y = mean))

b3 + geom_point(color="navy")+
  geom_smooth(method = "lm", color = "navy", fill = "lightgray")+
  stat_cor(method = "pearson",label.x = 0, label.y = 4.5) 

