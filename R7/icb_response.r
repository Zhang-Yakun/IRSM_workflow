setwd("D:/workplace/geneset/FIG6/")

library(dplyr)
library(data.table)
library(stringr)
library(tidyverse)

# 药物反应数据预处理 ---------------------------------------------------------------

gse1<-fread("GSE93157_info.txt",stringsAsFactors=FALSE)%>%as.data.frame()

gse1[7,]%>%as.character()%>%table()
#肺鳞癌=13，其他肺癌=22，黑色素瘤=25，头颈=5

gse2<-gse1[,gse1[7,]=="LUNG NON-SQUAMOUS CANCER"]
gse3<-gse1[,gse1[7,]=="SQUAMOUS LUNG CANCER"]

gse<-cbind(gse2,gse3)
colnames(gse)<-gse[1,]

gse[20,]%>%as.character()%>%table()
#反应数据:CR=1,PD=14,PR=8,SD=12
#CR, complete response; PR, partial response; SD, stable disease; PD, progression disease

gse[19,]%>%as.character()%>%table()
#单抗药：NIVOLUMAB=18，PEMBROLIZUMAB =17

mydata<-gse[c(7,19,20,9,18),]
mydata[6,]<-ifelse(mydata[3,]=="best.resp PD","PD","NPD")

data_n<-mydata[,mydata[2,]=="drug NIVOLUMAB"]
data_p<-mydata[,mydata[2,]=="drug PEMBROLIZUMAB"]


# 表达谱 ---------------------------------------------------------------------

riskgene<-fread("lasso_gene_index.CSV")%>%as.data.frame()

exp<-fread("GSE93157_series_matrix.txt")%>%as.data.frame()

rownames(exp)<-exp$ID_REF
exp<-exp[,-1]

exp_risk<-exp[riskgene$geneids,]

exp_risk2<-exp_risk[,colnames(mydata)]

group_response<-data.frame(sample=colnames(mydata),
                           drug=as.character(mydata[2,]),
                           response=as.character(mydata[6,]),
                           age=as.numeric(str_remove(as.character(mydata[4,]),pattern = "[a-z]* ")),
                           sex=as.character(mydata[5,]))
identical(colnames(exp_risk2),group_response$sample)

#####计算风险得分

gene<-riskgene$geneids
coef<-riskgene$index.min

risk.score<-c()
for(i in 1:ncol(exp_risk2)){
  risk.score[i]<-sum(coef*exp_risk2[,i])
}

medi<-median(risk.score)

rs.df<-data.frame(sample=colnames(exp_risk2),risk.score=risk.score,risk.group=ifelse(risk.score>medi,"High_risk","Low_risk"))

res_df<-inner_join(rs.df,group_response)

# heatmap --------------------------------------------------------------

res_df<-arrange(res_df,risk.group,response)

res_df$Age<-ifelse(res_df$age>=60,">=60","<60")

anno_col = data.frame(
  Group = as.factor(res_df$risk.group),
  Response = as.factor(res_df$response),
  Drug = as.factor(res_df$drug),
  Age = as.factor(res_df$Age),
  sex = as.factor(res_df$sex)
)

rownames(anno_col) = res_df$sample

exp_icb<-exp[c("PDCD1"),res_df$sample]

heat_input<-rbind(exp_risk2[,res_df$sample],exp_icb)

library(pheatmap)
pheatmap(heat_input,cellwidth = 6, cellheight = 14, fontsize = 8,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("navy", "white", "firebrick2"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = anno_col,
         treeheight_col = "0",#不画树
         border_color = "NA")


# violinplot -----------------------------------------------------------------

library(ggplot2)

e<-t(heat_input)%>%as.data.frame()
e$sample<-rownames(e)

vio_input<-inner_join(res_df,e)
vio_input$response<-as.factor(vio_input$response)
vio_input$risk.group<-as.factor(vio_input$risk.group)

######
#"CD40LG" "CX3CR1" "IL7"    "TLR3"   "PDCD1"  "TLR2"

ggplot(vio_input, aes(x=response, y=CD40LG,fill=response)) + 
    geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
    #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
    geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
    scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                   size=16),
          legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                    size=18),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+  #不显示网格线
    ylab("CD40LG")+xlab("") #设置x轴和y轴的标题


###
t.test(vio_input$CD40LG~vio_input$response)[[3]]
t.test(vio_input$CX3CR1~vio_input$response)[[3]]
t.test(vio_input$IL7~vio_input$response)[[3]]
t.test(vio_input$TLR3~vio_input$response)[[3]]
t.test(vio_input$PDCD1~vio_input$response)[[3]]
t.test(vio_input$TLR2~vio_input$response)[[3]]


# xgboost -----------------------------------------------------------------
library(xgboost)
library(caret)
library(ggplot2)
library(lattice)
library(pROC)

####输入数据预处理

expr<-t(exp_risk2)%>%as.data.frame()
expr$sample=rownames(expr)
expr2<-inner_join(expr,group_response)

input<-expr2[,c(1:5,8)]
input$response<-as.factor(input$response)
#表达谱和标签PD,NPD

trainlist<-createDataPartition(input$response,p=0.8,list=F)
trainset<-input[trainlist,]
testset<-input[-trainlist,]
##分训练集和验证集

####### 训练集的数据预处理

library(Matrix)
# 将trainset的1-4列（自变量）转换为矩阵
traindata1 <- data.matrix(trainset[,c(1:5)])
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse = T)
# 将因变量转换为numeric类型，-1是为了从0开始计数
train_y <- as.numeric(trainset[,6])-1
# 将自变量和因变量拼接为list
traindata <- list(data=traindata2,label=train_y) 
# data是自变量，label是因变量

# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label) 

####### 测试集的数据预处理

# 将自变量转化为矩阵
testset1 <- data.matrix(testset[,c(1:5)])

# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
testset2 <- Matrix(testset1,sparse=T) 
# 将因变量转化为numeric
test_y <- as.numeric(testset[,6])-1
# 将自变量和因变量拼接为list
testset <- list(data=testset2,label=test_y)
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data=testset$data,label=testset$label)

### 建立模型
model_xgb <- xgboost(data=dtrain,booster='gbtree',max_depth=6,eta=0.5,objective='multi:softmax',num_class=3,nround=25)
#用测试集预测
pre <- predict(model_xgb,newdata=dtest)
#模型评估
library(caret)
xgb.cf <-caret::confusionMatrix(as.factor(pre),as.factor(test_y))
xgb.cf

#####
#建立模型
bst <- xgboost(data=traindata$data,label=traindata$label,max.depth = 2,eta = 0.5,nround = 2,objective="binary:logistic")
#进行预测
pred <- predict(bst,testset$data)

#做交叉验证的函数参数与训练函数基本一致，只需要在原有参数的基础上设置nfold
cv.res <- xgb.cv(data = traindata$data,label=traindata$label,max.depth=2,eta=1,nround=2,objective="binary:logistic",nfold=5)

#####ROC曲线和AUC值
#在测试集上预测
pred2=round(predict(bst,newdata = testset$data))
#输出混淆矩阵
table(testset$label,pred2,dnn = c("真实值","预测值"))

##
require(pROC)
xgboost_roc<-roc(testset$label,as.numeric(pred2))
#绘制ROC曲线和AUC值
plot(xgboost_roc,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),gird.col=c("green","red"),max.auc.polygon=TRUE,auc.polygon.col="skyblue",print.thres=TRUE,main='xgboost模型和ROC曲线')


# xgboost2 ----------------------------------------------------------------

library(groupdata2)

split_matrix_cv <- function(inputArr=NULL, partition1_size=0.8, cat_col="PFS_status", num_col="PFS_days", id_col="PID", times=5){
  # 两种方法
  # R package groupdata2
  # https://rdrr.io/github/LudvigOlsen/R-splitters/man/partition.html
  # R package caret
  # https://www.rdocumentation.org/packages/caret/versions/6.0-86/topics/createDataPartition
  library(groupdata2)
  res_list <- list()
  for(i in 1:times){
    res_list_tmp <- partition(data=inputArr, p=partition1_size, cat_col=cat_col, num_col=num_col, id_col=id_col)
    names(res_list_tmp) <- c("train", "validate")
    res_list[[i]] <- res_list_tmp
  }
  return(res_list)
}

machine_learning <- function(inputArr=NULL, variable_cols=NULL, response_col=NULL, method=c("randomForest", "xgboost"), seed_num=1234, eval_metric = "merror", objective = "multi:softmax"){
  machine_learning_obj <- NA
  # 训练数据提取，提取因变量factor
  response_variable <- as.factor(as.character(inputArr[, response_col]))
  # 训练数据提取，提取自变量矩阵
  train_arr <- inputArr[, variable_cols]
  # 随机森林
  if(method[1]=="randomForest"){
    library(randomForest)
    machine_learning_obj <- randomForest(y=response_variable, x=train_arr)
  }
  if(method[1]=="xgboost"){
    library(xgboost)
    set.seed(seed_num)
    machine_learning_obj <- xgboost(data = as.matrix(inputArr[, variable_cols]), label = inputArr[, response_col], 
                                    eta = 0.3, max_depth = 15, nround=25, subsample = 0.5, colsample_bytree = 0.5,
                                    eval_metric = eval_metric, objective = objective, nthread = 24)
  }
  return(machine_learning_obj)
}


machine_learning_cv <- function(inputArr=NULL, variable_cols=NULL, response_col=NULL, id_col="PID", partition1_size=0.8, fold=5, 
                                method=c("randomForest"), eval_metric = "merror", objective = "multi:softmax", out_dir=NULL, out_PR=FALSE){
  # 构建交叉证实数据list 
  data_cv_list <- split_matrix_cv(inputArr=inputArr, partition1_size=partition1_size, cat_col=response_col, num_col=NULL, id_col=id_col, times=fold)
  # 显示数据分割情况
  print(lapply(data_cv_list, function(x){
    (sapply(x, dim))
  }))
  # split_matrix_cv(inputArr=inputArr, fold=5)
  # 构建返回模型对象list
  model_obj_list <- list()
  # 构建AUC list
  res_df_list <- list()
  for(i in 1:length(data_cv_list)){
    train_data <- as.data.frame(data_cv_list[[i]]$train, check.names=FALSE)
    validate_data <- as.data.frame(data_cv_list[[i]]$validate, check.names=FALSE)
    # if(method[1]=="randomForest"){
    # library(randomForest) 
    # 构建机器学习模型
    model_obj <- machine_learning(inputArr=train_data, variable_cols=variable_cols, response_col=response_col, method=method, eval_metric = eval_metric, objective = objective)
    # 对验证集进行预测
    res_pred <- predict(model_obj, as.matrix(validate_data[, variable_cols]), type="prob")
    # 提取验证集样本标签
    test_label <- as.character(validate_data[, response_col])
    # 计算AUC
    outprefix <- NULL
    # 是否输出AUC图片
    if(!is.null(out_dir)){
      outprefix <- paste0(out_dir, "model_", i)
      checkDir(dirname(outprefix))
    }
    if(method=="randomForest"){
      res_df <- summary_AUC(test_label=test_label, predict_value=res_pred[, 2], outprefix=outprefix, out_PR=out_PR)
    }
    if(method=="xgboost"){
      res_df <- summary_AUC(test_label=test_label, predict_value=res_pred, outprefix=outprefix, out_PR=out_PR)
    }   
    res_df <- data.frame(Model=paste0("model_", i), res_df)
    # 构建返回值
    model_obj_list[[i]] <- model_obj
    res_df_list[[i]] <- res_df
    # }
  }
  # 合并AUC list为矩阵
  res_AUC_df <- Reduce("rbind", res_df_list)
  # 构建返回list
  return(list(model_obj_list=model_obj_list, res_AUC_df=res_AUC_df))
}

input2<-expr2[,c(1:4,5,7)]

input2$response<-ifelse(input2$response=="PD","0","1")
input2$response<-as.factor(input2$response)

input2$sample<-as.factor(input2$sample)
var_cols<-colnames(input2)[1:4]

res_xgboost_cv <- machine_learning_cv(inputArr=input2, variable_cols=var_cols, response_col="response", 
                                      id_col="sample",partition1_size=0.75, fold=5, method="xgboost", eval_metric = "rmse", objective = "binary:logistic", out_dir=outprefix_tmp)
print(res_xgboost_cv$res_AUC_df)
save(res_xgboost_cv, file=paste0(res_home, gset_ICB$GEO_ID[1] ,"res_xgboost_cv.Rdata"))
