library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library()
setwd("E:/Rwork/lung/R5/")

g<-read.csv("lasso_gene_index.csv")
g$geneids
#12个risk gene


# 风险基因富集分析 ----------------------------------------------------------------

KEGG <- read.csv("E:/Rwork/lung/HSA_KEGG.csv")
TERM2GENE <- KEGG[, c("KEGGID", "ENTREZID")]
TERM2NAME <- KEGG[, c("KEGGID", "DESCRIPTION")]

ID <- bitr(g$geneids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
yy <- enricher(ID$ENTREZID, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
head(yy)
yy1 <- setReadable(yy, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
yy2 <- yy1[yy1$pvalue < 0.05, asis = T]
head(yy2)

col=c("#D1E8ED","#B8D9DE","#9ED1D3","#97D3D3","#8AC9C7","#81CAC6","#86C6C4","#71ADA8","#66A7A4","#589A96")

kegg<-yy2@result%>%arrange(p.adjust)

library(ggplot2)
ggplot(data=kegg[1:10,],aes(x=reorder(as.factor(Description),-log10(pvalue)),y=-log10(pvalue),fill=pvalue))+
  geom_bar(stat = "identity")+
  scale_fill_gradient(low="#589A96", high="#9ED1D3")+
  coord_flip()+
  labs(title="KEGG pathways",y="Enrichment Score",x="")


# 上调基因和风险基因的交集 ------------------------------------------------------------

up<-read.table("E:/Rwork/lung/老年组上调基因11.txt")
#11个上调基因

intersect(up$V1,g$geneids)
