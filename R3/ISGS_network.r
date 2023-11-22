setwd("D:/workplace/geneset/")

######90 gene
library(dplyr)
library(data.table)
g<-fread("immunosense_gene.txt",header=F)
g<-g$V1


# 基因circle图 ---------------------------------------------------------------

#####circlize
devtools::install_github("tidyverse/tidyverse")
BiocManager::install("ComplexHeatmap")
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

#####RCircos

library(RCircos)

data(UCSC.HG38.Human.CytoBandIdeogram)
head(UCSC.HG38.Human.CytoBandIdeogram)
#第一列 染色体编号；第二列 染色体片段起始位点；第三列 染色体片段结束位点；第四列 染色体片段编号；第五列 染色体片段颜色
# 这个数据是RCicos内置的人类染色体信息，第四列和第五列信息用于展示染色体的核型。可以省略

# 设置染色体数据
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram

RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=10, tracks.outside=0 )
# chr.exclude=NULL;  设置不显示的染色体，如 c(1,3)          
# tracks.inside=10;  设置内部环形个数
# tracks.outside=0;   设置外部环形个数 

# 绘制染色体图形，默认方法显示染色体名称。
RCircos.Set.Plot.Area()     
RCircos.Chromosome.Ideogram.Plot() 

####UCSC HG38 GTF

gtf<-fread("hg38_gene.bed")%>%as.data.frame()
gtf<-gtf[,c(1:3,6)]
colnames(gtf)<-c("Chromosome","chromStart","chromEnd","Gene")

RCircos.Gene<-subset(gtf,gtf$Gene%in%g)
rownames(RCircos.Gene)<-RCircos.Gene$Gene


# 指定内容在内侧的环形还是外侧的环形生成
side <- "in";

# 指定内容在第几个环形生成
track.num <- 1;
# 绘图
RCircos.Gene.Connector.Plot(RCircos.Gene, track.num, side);
# 在染色体上添加基因名称， 指定内容在第几个环形生成
name.col <- 4;
track.num <- 2;
RCircos.Gene.Name.Plot(RCircos.Gene, name.col,track.num, side);

####基因间的联系曲线

pp<-fread("string_interactions.tsv")%>%as.data.frame()

library(stringr)
pp$`#node1`<-str_replace(pp$`#node1`,"GNB2L1","RACK1")
pp$`node2`<-str_replace(pp$`node2`,"GNB2L1","RACK1")

pp<-subset(pp,pp$experimentally_determined_interaction!=0.000)

node<-pp[,1:2]
node1<-node$`#node1`%>%as.data.frame()
colnames(node1)<-"Gene"
node11<-inner_join(node1,gtf)
node11<-node11[,-1]

node2<-node$`node2`%>%as.data.frame()
colnames(node2)<-"Gene"
node22<-inner_join(node2,gtf)
node22<-node22[,-1]

link<-cbind(node11,node22)
colnames(link)[4:6]<-c("Chromosome.1","chromStart.1","chromEnd.1")

# 指定图形在第11个环形生成
track.num <- 6;
# 绘图
RCircos.Link.Plot(link, track.num, by.chromosome=TRUE);

##### GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

myID <- bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db,drop = TRUE)
ego <- enrichGO(
  gene=myID$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "ALL", # 或MF或CC
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 1
)
dim(ego)
go_res<-ego@result
write.csv(go_res,"GO-result.csv")

# GO柱状图 -------------------------------------------------------------------

library(stringr)
df <- go_res %>% dplyr::select(Description, ONTOLOGY, pvalue)
df$ONTOLOGY <- str_c("GO-", df$ONTOLOGY)

names(df) <- c("Name", "Pathway", "Score")
df <- df %>% group_by(Pathway) %>% dplyr::slice(1:10)
df$Score <- -log10(df$Score)

library(ggpubr)
library(ggthemes)
library(ggplot2)

ggbarplot(df, x = "Name", y="Score", fill = "Pathway", rotate = T, 
          palette = c("#ADC9D7", "#E39E3E", "#C7533B", "#509A80"), 
          xlab = "", ylab = "-log10(p-value)",
          sort.val = "asc")+
  theme_classic ()+
  #scale_y_continuous(limits = c(0, 20))+
  labs(fill ="Terms")


# 功能网络图 -------------------------------------------------------------------

library(enrichplot)
ego2<-pairwise_termsim(ego)
emapplot(ego2)

# 网络度展示 -------------------------------------------------------------------

d<-fread("degree.txt")%>%as.data.frame()
dd<-d[order(-d[,2]),]

library(ggplot2)
library(ggpubr)

dd$group<-ifelse(dd$Degree>30,"hub","no-hub")

ggdotchart(dd,x="name",y="Degree",
           col="group",
           palette = c("#FC4E07","#00AFBB"),
           sorting = "ascending",
           add="segments",
           ggtheme=theme_pubr())

h<-subset(dd,dd$group=="hub")

g11<-read.table("老年组上调基因11.txt")

setdiff(g11$V1,h$name)
