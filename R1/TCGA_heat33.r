setwd("D:/workplace/geneset/FIG1/")

library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
#brewer.pal(8,"Purples")
#brewer.pal(11,"PiYG")
#brewer.pal(11,"PRGn")
NES<-fread("TCGA_gsea_res.csv")%>%as.data.frame()

nes<-arrange(NES,desc(NES))

genum<-fread("TCGA_gene_number.csv")%>%as.data.frame()
rownames(genum)<-genum$type
gen<-genum[nes$V1,]

identical(gen$type,nes$V1)

#####

df1<-data.frame(cancer=nes$V1,
                A=rep("NES",length(nes$V1)),
                S=nes$NES
                )
df1$cancer <- factor(df1$cancer,levels=rev(df1$cancer))
#
df2<-data.frame(cancer=nes$V1,
                A=rep("pvalue",length(nes$V1)),
                S=nes$p.adjust
                )
df2$cancer <- factor(df2$cancer,levels=rev(df2$cancer))
#
df3<-data.frame(cancer=gen$type,
                A=rep("up",length(gen$type)),
                S=gen$up_regulate
)
df3$cancer <- factor(df3$cancer,levels=rev(df3$cancer))
#
df4<-data.frame(cancer=gen$type,
                A=rep("common",length(gen$type)),
                S=gen$common
)
df4$cancer <- factor(df4$cancer,levels=rev(df4$cancer))

#####NES

ggplot(data=df1, aes(A, cancer)) + 
  geom_tile(aes(fill = S), colour = "white",size=0.5)+
  scale_fill_gradient2(low="#4393C3",mid = "white",high = "#D6604D")+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size=10),# 调整x轴文字
        axis.text.y = element_text(size = 10)) +#调整y轴文字
  labs(fill ="NES")

#####pvalue
ggplot(data=df2, aes(A, cancer)) + 
  geom_tile(aes(fill = S), colour = "white",size=0.5)+
  scale_fill_gradient2(low = "#6A51A3",mid = "#9E9AC8",high = "#EFEDF5")+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size=10),# 调整x轴文字
        axis.text.y = element_text(size = 10)) +#调整y轴文字
  labs(fill ="pvalue")

#####up_gene

ggplot(data=df3,aes(A,cancer)) + 
  geom_point(aes(color = S,size=as.factor(S)))+
  scale_color_gradient2(low = "#D9F0D3",mid = "#A6DBA0",high = "#00441B")+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size = 10),# 调整x轴文字
        axis.text.y = element_text(size = 10))+#调整y轴文字
  #调整legend
  labs(color ="old specific gene")

ggplot(data=df4,aes(A,cancer)) + 
  geom_point(aes(color = S,size=as.factor(S)))+
  scale_color_gradient2(low = "#FDE0EF",mid = "#F1B6DA",high = "#8E0152")+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size = 10),# 调整x轴文字
        axis.text.y = element_text(size = 10))+#调整y轴文字
  #调整legend
  labs(color ="common gene")

# ladar -------------------------------------------------------------------

library(ggradar)
library(dplyr)
library(scales)
library(tibble)
library(data.table)
library(reshape2)
library(ggplot2)

setwd("D:/workplace/geneset/FIG1/gene_radar/TCGA")
glist<-list.files("D:/workplace/geneset/FIG1/gene_radar/TCGA")

####tcga

for(i in 1:3){
  
  g1<-fread(glist[i])%>%as.data.frame()
  ra1<-g1%>%reshape2::acast(gene~cancer)%>%as.data.frame()
  ra1$group<-rownames(ra1)
  
  ra1<-ra1%>%select(group,everything())
  ####
  ggradar(ra1, grid.min = -3,
          grid.mid = 0, grid.max = 3,
          values.radar = c("-3", "0", "3"),
          gridline.min.colour = "grey",
          gridline.mid.colour = "blue", gridline.max.colour = "orange",
          axis.label.size = 5, axis.line.colour = "grey",
          plot.title = paste("log2(Fold Change) of", rownames(ra1), "in TCGA",sep = " "), legend.text.size = 14, legend.position = "left",
          background.circle.colour = "white",
          background.circle.transparency = 0.1,
  )
  ggsave(paste(rownames(ra1),".pdf",sep=""), width = 7.33, height = 4.54, units = "in")
}

#################

setwd("D:/workplace/geneset/FIG1/gene_radar/GTEX")
glist<-list.files("D:/workplace/geneset/FIG1/gene_radar/GTEX")

####GTEX
  i=3 
  g1<-fread(glist[i])%>%as.data.frame()
  ra1<-g1%>%reshape2::acast(gene~tissue)%>%as.data.frame()
  ra1$group<-rownames(ra1)
  
  ra1<-ra1%>%select(group,everything())
  ####
  ggradar(ra1, grid.min = -3,
          grid.mid = 0, grid.max = 15,
          values.radar = c("-3", "0", "15"),
          gridline.min.colour = "grey",
          gridline.mid.colour = "blue", gridline.max.colour = "orange",
          axis.label.size = 5, axis.line.colour = "grey",
          plot.title = paste("log2(Fold Change) of", rownames(ra1), "in GTEX",sep = " "), legend.text.size = 14, legend.position = "left",
          background.circle.colour = "white",
          background.circle.transparency = 0.1,
  )
  ggsave(paste(rownames(ra1),".pdf",sep=""), width = 7.33, height = 4.54, units = "in")


