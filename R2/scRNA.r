
library(Seurat)
library(data.table)
library(dplyr)
library(stringr)

setwd("E:/Rwork/lung/singlecell/")

rna.list = readRDS("rna_data.rds") 

head(rna.list@meta.data)
Idents(rna.list)
#原来细胞聚类的名称

new.cluster.ids<-paste(rna.list@meta.data$cell_type) 

names(new.cluster.ids)<-levels(rna.list)

rna.list<-RenameIdents(rna.list,new.cluster.ids)
Idents(rna.list)
head(rna.list@meta.data)

######UMAP

DimPlot(rna.list,reduction = "umap",pt.size = 0.1,group.by = "cell_type",label = TRUE)+NoLegend()
DimPlot(rna.list,reduction = "umap",label = TRUE,pt.size = 0.1)+NoLegend()

##### percent barplot

library(ggplot2)
library(reshape2)

#输入包含2个维度的数据。每一行是一簇柱子，每一列用不同的颜色图例区分

rt<-data.frame(sample=rna.list@meta.data$sample,celltype=rna.list@meta.data$cell_type)

df = reshape2::melt(table(rt))
# 把数据转换成ggplot常用的类型（长数据）

p = ggplot(df, aes(x=factor(sample,levels =unique(sample)),  # 转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
                   y=value, 
                   fill=factor(celltype,levels = unique(celltype)), 
))+
  labs(
    x="",   # 调整x轴名称
    y="",   # 调整y轴名称
    fill="" # 调整图例名称
  )
p + geom_bar(
  position="fill",
  stat="identity"
)+theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

######ssGSEA,这套数据不可用

library(escape)
library(GSEABase)
library(dittoSeq)
library(irGSEA)
library(UCell)

geneset<-clusterProfiler::read.gmt("E:/Rwork/lung/my_geneset.gmt")
myGS<-list(immunosenescence=geneset$gene)

#enrichScore = enrichIt(obj = rna.list, gene.sets = GS,method = "ssGSEA", min.size = NULL)

GS <- getGeneSets(library = "C7")
rna.list <- suppressWarnings(rna.list)
ES <- enrichIt(obj = rna.list, 
               gene.sets = GS, groups = 1000, cores = 2)

######addmoduleScore

geneset<-clusterProfiler::read.gmt("my_geneset.gmt")
myGS<-list(immunosenescence=geneset$gene)

myscore<-AddModuleScore(object = rna.list,myGS,name="immunosenescence")

head(myscore@meta.data)

DimPlot(myscore,reduction = "umap",pt.size = 0.1,group.by = "immunosenescence1")+NoLegend()
######

gene4<-c("CD40LG","CX3CR1","IL7","TLR3")

myGS2<-list(riskscore=gene4)

myscore2<-AddModuleScore(object = rna.list,myGS2,name="riskscore")

head(myscore2@meta.data)

#FeaturePlot(object=myscore,features = "immunosenescence1",reduction = "umap",cols = c("grey","red")) 

cell<-c("Treg","NK","Cytotoxic_CD8_T","CD4_T","B")

df2<-subset(myscore2@meta.data,
            myscore2@meta.data$cell_type=="B")
df1<-subset(myscore@meta.data,
            myscore@meta.data$cell_type=="B")

cor.test(df2$riskscore1,
         df1$immunosenescence1)

########直接比较两组样本的免疫衰老得分差异

score_df<-data.frame(sample=myscore@meta.data$sample,
                     cell=myscore@meta.data$cell_type,
                     cluster=myscore@meta.data$seurat_clusters,
                     score=myscore@meta.data$immunosenescence1)
rownames(score_df)<-rownames(myscore@meta.data)

age_df<-data.frame(sample=c("ERX2757104","ERX2757105","ERX2757106","ERX2757107",
                            "ERX2757108","ERX2757109","ERX2757111","ERX2757112",
                            "ERX2757113"),
                   age=c(65,65,65,60,60,60,55,55,55))
age_df$group<-ifelse(age_df$age>=60,"old","young")

old_score_df<-inner_join(score_df,age_df)
rownames(old_score_df)<-rownames(score_df)

#########boxplot
wilcox.test(score~as.factor(group),data = old_score_df,)

boxplot(score~as.factor(group),data = old_score_df,col=c("#F0CC8D","#C6DBC6"))
# "#6BAED6","#C6DBEF","#FCBBA1","#FB6A4A"
#old:"#F0CC8D",young:"#C6DBC6",c("#F0CC8D","#C6DBC6")

quantile(old_score_df$score)
old_score_df$scatter<-ifelse(old_score_df$score<0.02800432,"A",
                             ifelse(old_score_df$score<0.08538384,"B",
                                    ifelse(old_score_df$score<0.14213784,"C","D")))
old_score_df$scatter<-as.factor(old_score_df$scatter)
identical(rownames(myscore@meta.data),rownames(old_score_df))
myscore@meta.data$scatter<-old_score_df$scatter
DimPlot(myscore,reduction = "umap",pt.size = 0.1,group.by = "scatter",
        cols = c("#6BAED6","#FCBBA1","#FC9272","#FB6A4A") )

##########barplot

rt<-data.frame(sample=old_score_df$group,celltype=old_score_df$cell)

df = reshape2::melt(table(rt))
# 把数据转换成ggplot常用的类型（长数据）

library(ggsci)

p = ggplot(df, aes(x=factor(sample,levels =unique(sample)),  # 转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
                   y=value, 
                   fill=factor(celltype,levels = unique(celltype)), 
))+
  labs(
    x="",   # 调整x轴名称
    y="",   # 调整y轴名称
    fill="" # 调整图例名称
  )
p + geom_bar(
  position="fill",
  stat="identity"
)+theme_bw()

#pal_d3("category20")(17)
# scale_fill_igv(),scale_fill_d3("category20")

###########免疫衰老相关细胞筛选

#识别差异基因
scRNA.pos.markers<-FindAllMarkers(myscore,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)

write.csv(scRNA.pos.markers,file = "E:/Rwork/lung/singlecell/Allmarkers-pos.csv")
#该markers文件用于后边的细胞类型注释
head(scRNA.pos.markers)

Top_pos_genes<-scRNA.pos.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
#找出各个cluster中的TOP差异基因，按logFC值排序，n=2展示各个cluster的TOP2



# 山脊图-每类免疫细胞的score --------------------------------------------------------

ES<-myscore@meta.data

library(reshape2)
library(ggplot2)
library(ggridges)

ggplot(ES,aes(x=immunosenescence1,y=cell_type,fill=cell_type))+
  geom_density_ridges(alpha=0.5)+
  theme_ridges()+
  theme(legend.position = "none")

# upset图--得分高的细胞之间的上调基因 -------------------------------------------------------------------

###得分高的细胞类型

cell<-c("Treg","NK","Cytotoxic_CD8_T","CD4_T","B")

####upset输入数据
scRNA.pos.markers<-read.csv("Allmarkers-pos.csv")
t<-table(scRNA.pos.markers$cluster,scRNA.pos.markers$gene)

tt<-t(t)
UpSet_input<-tt[,cell]%>%as.matrix.data.frame()
UpSet_input<-UpSet_input%>%as.data.frame()

rownames(UpSet_input)<-rownames(tt)
colnames(UpSet_input)<-c("Treg","NK","Cytotoxic_CD8_T","CD4_T","B")

knitr::kable(head(UpSet_input))

####可视化

#devtools::install_github("hms-dbmi/UpSetR")
library(UpSetR)

upset(UpSet_input, sets = c("Treg","NK","Cytotoxic_CD8_T","CD4_T","B"),
      sets.bar.color = "#56B4E9",order.by = "freq",
      empty.intersections = "on")

# CD4_T和NK 共有上调基因的功能富集 ----------------------------------------------------

gene_NK <- subset(UpSet_input,UpSet_input$NK==1)
common_gene <- subset(gene_NK,gene_NK$CD4_T==1)
common_gene89<-rownames(common_gene)[-1]

####89个共有基因进行KEGG功能富集

library(clusterProfiler)
library(org.Hs.eg.db)

# KEGG
KEGG <- read.csv("HSA_KEGG.csv")
TERM2GENE <- KEGG[, c("KEGGID", "ENTREZID")]
TERM2NAME <- KEGG[, c("KEGGID", "DESCRIPTION")]

ID <- bitr(common_gene89, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
yy <- enricher(ID$ENTREZID, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
head(yy)
yy1 <- setReadable(yy, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(yy1)
yy2 <- yy1[yy1$p.adjust < 0.05, asis = T]
head(yy2) #14条显著富集的通路

write.csv(yy2,"NK_CD4T_KEGG.csv",row.names = F)

library(ggplot2)

yy3=data.frame(pvalue=yy2$pvalue,KEGG_pathway=yy2$Description,Count=yy2$Count)

yy3<-arrange(yy3,pvalue)
yy3$KEGG_pathway<-factor(yy3$KEGG_pathway,levels = rev(yy3$KEGG_pathway))
p<-ggplot(yy3,aes(pvalue,KEGG_pathway))
pbubble<-p+geom_point(aes(size=Count/2,color=(-1)*log10(pvalue)))
pbubble+scale_color_gradient(low="#C5CAE9FF",high="#273492FF")+labs(color=expression(-log[10](pvalue)),size="Count",
                                                          x="pvalue",y="KEGG pathway",title = "Pathway enrichment")+theme_bw()



# cellchat 年轻组 ---------------------------------------------------------------

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)

setwd("E:/Rwork/lung/singlecell/cellchat/")

#设置硬件参数，8线程
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))


# 数据读取

data.input <-  rna.list@assays$RNA # normalized data matrix

meta <- rna.list@meta.data
meta$group<-ifelse(meta$sample%in%c("ERX2757111","ERX2757112","ERX2757113"),"young","old")

unique(meta$group)
cell.use <- rownames(meta)[meta$group == "young"] 
# 按指定的变量提取细胞

data.input <- data.input[,cell.use]#取出对应细胞,也就是说，data的列名是meta的行名
data.input[1:5,1:5]

meta = meta[cell.use, ]#取出对应细胞的meta信息
unique(meta$cell_type)#看meta中储存的细胞注释信息，稍后用它作为分组依据

#####创建cellchat对象

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")

#创建celllchat对象，group.by指定通讯间的对象，用meta中的注释作为分组依据
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cell_type") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) #每种细胞的细胞数量


#######设置参考数据库

CellChatDB <- CellChatDB.human
#查看数据库的组成比例
showDatabaseCategory(CellChatDB)
# 查看数据库具体信息
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中，相当于设置好接下来要参考的数据库

##########表达数据预处理

cellchat <- subsetData(cellchat)#取出表达数据
cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因#
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- projectData(cellchat, PPI.human)#投影到PPI，储存上一步的结果到cellchat@LR$LRsig

#####推断受体配体信号网络

## 1.在配受体水平上计算细胞通讯
#计算细胞与细胞之间通信的概率
cellchat <- computeCommunProb(cellchat, raw.use = T)#默认计算方式为#type = "truncatedMean",
##去掉通讯数量很少的细胞，默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat <- filterCommunication(cellchat, min.cells = 10)

#将推断的结果提取出来
df.net <- subsetCommunication(cellchat)#将细胞通讯预测结果以数据框的形式取出
#install.packages("DT"),DT就是展示数据框的一种工具，可以调节展示的参数等
library(DT)
DT::datatable(df.net)
write.csv(df.net,'E:/Rwork/lung/singlecell/cellchat/young.df.net.csv')

## 2.在细胞通路水平上计算细胞间通讯，CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平的通信概率。

cellchat <- computeCommunProbPathway(cellchat)

## 3.计算聚合的细胞间通讯网络。我们可以通过计算链路的数量或汇总通信概率来计算聚合的细胞通信网络。用户还可以通过设置 sources.use和 targets.use来计算细胞组子集之间的聚合网络，使用circle plot显示任意两个细胞组之间的交互次数或总交互强度（权重）

cellchat <- aggregateNet(cellchat)#计算细胞-细胞聚合通信网络
groupSize <- as.numeric(table(cellchat@idents))
groupSize
#组图，1行3列
par(mfrow = c(1,3), xpd=TRUE)
#使用圆图显示任意两个细胞群之间的相互作用数量或总相互作用强度(权重)。颜色和source是一致的，圈圈的大小是每个细胞群细胞的数量
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
#总相互作用强度
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
#以luminal_macrophage为中心的相互作用强度
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'luminal_macrophage')


#检查每种细胞和其他细胞的互作情况
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}
dev.off()

saveRDS(cellchat, file = "cellchat_young_vs_old.rds")


# 老年组中-细胞通讯 ---------------------------------------------------------------

# 数据读取

data.input <-  rna.list@assays$RNA # normalized data matrix

meta <- rna.list@meta.data
meta$group<-ifelse(meta$sample%in%c("ERX2757111","ERX2757112","ERX2757113"),"young","old")

unique(meta$group)
cell.use <- rownames(meta)[meta$group == "old"] 
# 按指定的变量提取细胞

data.input <- data.input[,cell.use]#取出对应细胞,也就是说，data的列名是meta的行名
data.input[1:5,1:5]

meta = meta[cell.use, ]#取出对应细胞的meta信息
unique(meta$cell_type)#看meta中储存的细胞注释信息，稍后用它作为分组依据

#####创建cellchat对象

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")

#创建celllchat对象，group.by指定通讯间的对象，用meta中的注释作为分组依据
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cell_type") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) #每种细胞的细胞数量


#######设置参考数据库

CellChatDB <- CellChatDB.human
#查看数据库的组成比例
showDatabaseCategory(CellChatDB)
# 查看数据库具体信息
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中，相当于设置好接下来要参考的数据库

##########表达数据预处理

cellchat <- subsetData(cellchat)#取出表达数据
cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因#
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- projectData(cellchat, PPI.human)#投影到PPI，储存上一步的结果到cellchat@LR$LRsig

#####推断受体配体信号网络

## 1.在配受体水平上计算细胞通讯
#计算细胞与细胞之间通信的概率
cellchat <- computeCommunProb(cellchat, raw.use = T)#默认计算方式为#type = "truncatedMean",
##去掉通讯数量很少的细胞，默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat <- filterCommunication(cellchat, min.cells = 10)

#将推断的结果提取出来
df.net <- subsetCommunication(cellchat)#将细胞通讯预测结果以数据框的形式取出
#install.packages("DT"),DT就是展示数据框的一种工具，可以调节展示的参数等
library(DT)
DT::datatable(df.net)
write.csv(df.net,'E:/Rwork/lung/singlecell/cellchat/old.df.net.csv')

## 2.在细胞通路水平上计算细胞间通讯，CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平的通信概率。

cellchat <- computeCommunProbPathway(cellchat)

## 3.计算聚合的细胞间通讯网络。我们可以通过计算链路的数量或汇总通信概率来计算聚合的细胞通信网络。用户还可以通过设置 sources.use和 targets.use来计算细胞组子集之间的聚合网络，使用circle plot显示任意两个细胞组之间的交互次数或总交互强度（权重）

cellchat <- aggregateNet(cellchat)#计算细胞-细胞聚合通信网络
groupSize <- as.numeric(table(cellchat@idents))
groupSize
#组图，1行3列
par(mfrow = c(1,3), xpd=TRUE)
#使用圆图显示任意两个细胞群之间的相互作用数量或总相互作用强度(权重)。颜色和source是一致的，圈圈的大小是每个细胞群细胞的数量
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
#总相互作用强度
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
#以CD4_T为中心的相互作用强度
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'CD4_T')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'B')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'Cytotoxic_CD8_T')

#检查每种细胞和其他细胞的互作情况
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}
dev.off()

######可视化信号通路网络

pathways.show <- df.net$pathway_name#计算到的所有通路
#如果要选择其中一个信号的话，pathways.show <- c("TGFb")

#层次图可视化信号通路,左侧部分显示对某些感兴趣的细胞（即定义的vertex.receiver）的自分泌和旁分泌信号，右侧部分显示对数据集中剩余细胞的自分泌和旁分泌信号，在层次图中，实心圆和空心圆分别代表源和目标。

#用户应定义vertex.receiver，这是一个数值向量，代表的是感兴趣的细胞类型，比如(1,4)是指在Hierarchy plot中，左侧的感兴趣的细胞群的信号传导，这些感兴趣的细胞群是cellchat@idents里细胞的前四种，也就是df.net里target的前4种。
#CD4_T, B cell, Cytotoxic_CD8_T, Plasma_B
levels(cellchat@idents)

vertex.receiver = c(4,5,8,16)
netVisual_aggregate(cellchat, signaling =pathways.show[1],  
                    vertex.receiver = vertex.receiver,layout = 'hierarchy')
