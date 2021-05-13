setwd("/home/data/human/ESCC/Linna-ESCC/Fig/source_data/")
library(scales)
library(ggplot2)

####CD45- TSNE####
df1 <- readRDS("/home/data/human/ESCC/Linna-ESCC/Fig/Cellinfo_20190924/Cell.info/CD45-_addTEC5_info.rds")
df_EPI=df1[df1$celltype=="Epithelial",]
df_FIBRO=df1[df1$celltype=="Fibroblast",]
df_Endo=df1[df1$celltype=="Endothelial",]
df_FRC=df1[df1$celltype=="FRC",]
df_peri=df1[df1$celltype=='Pericytes',]

col3<-c("#f57665",'#1279A2',"#CAA57D","#f1c550","#0b8457")

ggplot(df_FRC, aes(tSNE1, tSNE2))+
  geom_point(size=0.3,alpha=0.6,color="#0b8457")+
  geom_point(data=df_EPI, aes(tSNE1, tSNE2), color="#f57665", size=0.3,alpha=0.6)+
  geom_point(data=df_Endo, aes(tSNE1, tSNE2), color="#f1c550",size=0.3,alpha=0.6)+
  geom_point(data=df_FIBRO, aes(tSNE1, tSNE2), color="#1279A2", size=0.3,alpha=0.6)+
  geom_point(data=df_peri, aes(tSNE1, tSNE2), color="#CAA57D", size=0.3,alpha=0.6)+
  guides(color=guide_legend(title=NULL,override.aes = list(size = 5)))+
  theme_linedraw()+theme(panel.grid =element_blank())+
  theme(legend.position = 'none' )
#ggsave("CD45-_TSNE.png",res1, height = 6, width = 6)

ggplot(df1, aes(tSNE1, tSNE2,color=celltype))+
  geom_point(size=0.3)+ scale_color_manual(values=col3)+
  guides(color=guide_legend(title=NULL,override.aes = list(size = 5)))+
  theme_linedraw()+theme(panel.grid =element_blank())+
  theme(legend.position="bottom")
#ggsave("CD45-_TSNE_legend.png",res2, height = 6, width = 6)

df_t=df1[df1$tissue=='Tumor',]
df_n=df1[df1$tissue=='Adj. normal',]

ggplot(df_t, aes(tSNE1, tSNE2))+
  geom_point(size=0.3,alpha=0.8,color="#C6DBEF")+
  geom_point(data=df_n, aes(tSNE1, tSNE2), color="#08306B", size=0.3,alpha=0.8)+
  guides(color=guide_legend(title=NULL,override.aes = list(size = 5)))+
  theme_linedraw()+theme(panel.grid =element_blank())+
  theme(legend.position = 'none' )


ggplot(df1, aes(tSNE1, tSNE2, color=tissue))+
  geom_point(size=0.3)+ scale_color_manual(values=rev(c("#C6DBEF","#08306B")))+
  guides(color=guide_legend(title=NULL,override.aes = list(size = 5)))+
  theme_linedraw()+theme(panel.grid =element_blank())+
  theme(legend.position="bottom")

#ggsave("CD45-_Tissue.png",res3, height = 6, width = 6)
#ggsave("CD45-_Tissue_legend.png",res4, height = 6, width = 6)
edit(colnames(df1))
result=df1[,c("cell", "tissue", "sample", "celltype",   "tSNE1", "tSNE2")]
write.table(result,"Fig1b.txt",quote = F,sep="\t",row.names = F)
