setwd("D:/workplace/geneset/FIG1/")

NES<-fread("GTEX_gsea_res.csv")%>%as.data.frame()

nes<-arrange(NES,desc(NES))

genum<-fread("GTEX_gene_number.csv")%>%as.data.frame()
rownames(genum)<-genum$type
gen<-genum[nes$V1,]

identical(gen$type,nes$V1)
####

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
####
df5<-data.frame(cancer=df1$cancer,
                NES=df1$S,
                p=df2$S,
                up=df3$S,
                com=df4$S)
write.csv(df5,"GTEX_29_all.csv")
