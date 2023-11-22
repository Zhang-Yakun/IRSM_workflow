

# 提取Msigdb的c7基因集 - immunologic signatures----------------------------------------------------------

library(msigdbr)

m_df = msigdbr(species = "Homo sapiens", category = "C7")

#利用split函数把得到的基因集分组建成一个list

m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
str(m_list)


# 免疫基因集的GSVA --------------------------------------------------------------

gsva_res_c7<- gsva(as.matrix(tpm_exp), m_list)
head(gsva_res_c7)  #GSVA result

gsva_res_7<-t(gsva_res_c7)%>%as.data.frame()

