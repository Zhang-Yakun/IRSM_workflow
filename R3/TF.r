setwd("E:/Rwork/lung/TF/")

###
immune_gene90<-read.table("immunosense_gene.txt",header=F) #90个免疫衰老基因列表
immune_gene90<-immune_gene90$V1%>%unique()%>%as.character()

# 三个转录因子数据库 ---------------------------------------------------------------

######TRRUST数据
TF1<-fread("trrust_rawdata.human.tsv")%>%as.data.frame()

myTF1<-subset(TF1,TF1$V2 %in% immune_gene90)
myTF1<-myTF1[,1:2]
colnames(myTF1)<-c("TF","target")

#######RegNetwork数据

TF2<-fread("human.source")%>%as.data.frame()

myTF2<-subset(TF2,TF2$V3 %in% immune_gene90)
myTF2<-myTF2[,c(1,3)]
colnames(myTF2)<-c("TF","target")

########hTFtarget数据

TF3<-fread("hTFtarget_TF-Target.txt")%>%as.data.frame()

myTF3<-subset(TF3,TF3$target%in% immune_gene90)
myTF3<-myTF3[,1:2]

#######合并三种数据

all_TF<-rbind(myTF1,myTF2)
all_TF<-rbind(all_TF,myTF3)
all_TF<-distinct(all_TF)

nrDEG<-unique(all_TF$target)

#88个被TF调控的免疫衰老基因

# 计算TF和靶基因的spearman相关性 ----------------------------------------------------

#读入表达数据
exp<-fread("exp_261.csv")%>%as.data.frame()
exp2<-exp[,5:dim(exp)[2]]
row.names(exp2)<-exp$V1

expt<-t(exp2)%>%as.data.frame()
TF_cor_df<-data.frame(TF=0,target=0,coefficient=0)

exptt<-expt[!apply(expt,1,sum)==0,]
nrDEG<-intersect(nrDEG,rownames(exptt))
#82个有表达值的被转录因子调控的免疫衰老基因

for(i in 1:length(nrDEG)){
  gene<-nrDEG[i]
  ge_TF<-subset(all_TF,all_TF$target==gene)
  
  TF_cor<-intersect(ge_TF$TF,rownames(exptt))
  for(j in 1:length(TF_cor)){
    TF<-TF_cor[j]
    x=exptt[gene,]%>%as.numeric()
    y=exptt[TF,]%>%as.numeric()
    coef<-cor(x,y,method = "spearman")
    
    if(abs(coef)>0.6){
      
      TF_cor_res<-data.frame(TF=TF,target=gene,coefficient=coef)
      TF_cor_df<-rbind(TF_cor_df,TF_cor_res)
    }
  }
  
}
TF_cor_df<-TF_cor_df[-1,]
TF_cor_df$TF%>%unique()%>%length()
TF_cor_df$target%>%unique()%>%length()
##得到35个TF，39个target的调控网络


write.csv(TF_cor_df,"TF-target-cor.csv",row.names = F)

