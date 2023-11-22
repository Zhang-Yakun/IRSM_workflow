
setwd("D:/workplace/geneset/geo/")

library(dplyr)
library(data.table)
library(ggplot2)

######风险评分绘图

### 读入风险评分数据

rs.df<-read.csv("RiskScore_GSE26939.csv")

### 2.风险评分绘图
library(ggplot2)
## 样品按score值排序
rs.df <- rs.df[order(rs.df$risk.score),]
rs.df2<-subset(rs.df,!is.na(rs.df$time))

p1 <- ggplot(data = rs.df2) +
  geom_point(aes(x=seq(1:length(rs.df2$X)), y=risk.score, color=risk.group)) + 
  scale_color_manual(values = c("#39478A","#DF1815")) + # 自定义颜色映射
  ggtitle("Risk Score") + 
  theme(plot.title = element_text(hjust = 0.5)) + # title居中
  xlab("") + 
  ylab("Risk Score") +
  geom_vline(aes(xintercept=(length(rs.df2$X))/2),colour = "#BB0000",linetype = "dashed")

### 3.生存状态散点图
p2 <- ggplot(data = rs.df2) +
  geom_point(aes(x=seq(1:length(rs.df2$X)), y=time, color=as.factor(OS))) + 
  ggtitle("Survival status") + 
  scale_color_manual(values = c("#39478A","#DF1815")) +
  theme(plot.title = element_text(hjust = 0.5)) + # title居中
  xlab("patients") + 
  ylab("Survival time (years)")

### 4.组成图
library(cowplot)
plot_grid(p2,p1,ncol = 1, align = "h",labels = c("A","B"))


# 表达热图 --------------------------------------------------------------------

exp<-read.csv("D:/workplace/geneset/geo/GSE26939/exp_GSE26939.csv")

rownames(exp)<-exp$X
exp<-exp[,-1]

exp2<-exp[,rs.df2$X]

####
anno_col = data.frame(
  Group = as.factor(rs.df2$risk.group)
)

rownames(anno_col) = rs.df2$X

library(pheatmap)
pheatmap(exp2,cellwidth = 2, cellheight = 18, fontsize = 8,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("navy", "white", "firebrick2"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = anno_col,
         treeheight_col = "0",#不画树
         border_color = "NA")
