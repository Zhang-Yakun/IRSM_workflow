
setwd("E:/Rwork/lung/R5/")

library(data.table)
library(dplyr)

df1<-read.csv("gsva_risk_score.csv") #两个得分

df2<-fread('exp_261.csv')%>%as.data.frame() #基因表达
rownames(df2)<-df2$V1

checkpoint<-c("PDCD1","CTLA4","TIGIT","LAG3")
df2<-df2[,checkpoint]
df2$sample<-rownames(df2)

check_df<-inner_join(df1,df2)


# 热图 ----------------------------------------------------------------------


# 相关性图 --------------------------------------------------------------------

cor(check_df[,7:12],method = 'spearman')

library(Hmisc) 
df_rcorr<-rcorr(as.matrix(check_df[,7:12])) 
#半角圆圈图
library(corrplot)
corrplot(df_rcorr$r,#数据
         type="lower",#可选择展示方式，"full", "lower", "upper"
         tl.col ="red",#文本颜色
         tl.srt = 45#标签旋转
)

corrplot(corr =df_rcorr$r, p.mat = env.p,method = "circle",type = "upper")
corrplot(corr =df_rcorr$r, type="lower",add=TRUE,method="number")


corrplot(
  # 相关系数矩阵
  corr = df_rcorr$r, 
  order = 'AOE',
  type = 'lower', 
  tl.pos = 'd'
  )


corrplot(
  corr = df_rcorr$r, 
  add = TRUE, 
  type = 'upper',  # 指定展示的方式，可以是完全的、下三角或上三角
  method = 'number',# 指定可视化的方法，可以是圆形、方形、椭圆形、数值、阴影、颜色或饼图形
  order = 'AOE', # 指定相关系数排序的方法，可以是原始顺序(original)、特征向量角序(AOE)、第一主成分顺序(FPC)、层次聚类顺序(hclust)和字母顺序，一般”AOE”排序结果都比”FPC”要好
  diag = FALSE,  #是否展示对角线上的结果，默认为TRUE
  tl.pos = 'n', # 指定文本标签(变量名称)的位置，当type=full时，默认标签位置在左边和顶部(lt)，当type=lower时，默认标签在左边和对角线(ld)，当type=upper时，默认标签在顶部和对角线，d表示对角线，n表示不添加文本标签
  cl.pos = 'n')# 图例（颜色）位置，当type=upper或full时，图例在右表(r)，当type=lower时，图例在底部，不需要图例时，只需指定该参数为n
