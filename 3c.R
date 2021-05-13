library("stringr")
library(reshape2)
library(limma)
library(GSVA)
library(pheatmap)
setwd("/home/data/human/ESCC/Linna-ESCC/Fig/Figure_20191004/Fig3/")
info <- readRDS("CellInfo_EpiProgram_0.7.rds")
gsva <- readRDS("/home/data/human/ESCC/Linna-ESCC/EPI_rmNormal_rmLess100/EPI_rmNormal_rmLess100_PC20_GSVA_MarkerGenes.rds")
gsva[1:5,1:5]


#######metainfo######
metainfo=readRDS("/home/data/human/ESCC/Linna-ESCC/EPI_rmNormal_rmLess100/Fig/EPI_rmNormal_rmLess100Cellinfo.rds")
metainfo$cluster=gsub("C","E",metainfo$cluster)
metainfo$cluster=as.factor(metainfo$cluster)

col3<-c("#DEEAB2","#64B473","#2D553C","#A1D8C9","#487C76","#7AAF93","#D0E4C6",
        "#F3746C","#BB4B95","#F8BFAF","#F7DDD3","#66CDF6","#598198","#D5E7F7",
        "#69B3CE","#D6D5EB","#7B8CBC","#7674AE", "#E3CEE4",
        "#AFB2B6","#C9BDB2","#F5C785","#ECAECF","#E9A943","#CAA57D","#A79388",
        "#EACF68","#F6F3A7","#C45337","#86382B","#EADCE4",
        "#EE5276","#9E6CAB","#74507B")
names(col3)=c("E01", "E02", "E11", "E12", "E13", "E14", "E15", "E16", "E17", 
              "E18", "E09", "E20", "E03", "E21", "E22", "E23", "E24", "E25", 
              "E26", "E27", "E28", "E29", "E30", "E04", "E05", "E06", "E07", 
              "E08", "E10", "E19")

ggplot(metainfo, aes(tSNE1, tSNE2, color=cluster))+
  geom_point(size=0.3)+scale_color_manual(values = col3)+
  guides(color=guide_legend(title=NULL,override.aes = list(size = 5)))+
  theme_linedraw()+theme(panel.grid =element_blank())+
  theme(legend.position = 'right' )+
  theme(aspect.ratio = 1)

setwd("/home/data/human/ESCC/Linna-ESCC/Fig/Figure_20191128/Fig3/")
######canculate P and T #######
info<-info[colnames(gsva),]
index=c("Innate","Stress","AP", "Cycling","Epi1","Epi2","Mes","Oxd")
est=gsva

result<-NULL
  
for (i in index){
  change=info[,i]
  design<-ifelse(change=="TRUE","A",'other')
  design <- model.matrix(~ 0+factor(design))
  head(design)
  colnames(design) <- c('A','ohter')
  head(design)
  contrast.matrix <- makeContrasts(A-ohter, levels=design)

fit <- lmFit(est, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
sigPathways <- topTable(fit2,  number=Inf,  adjust="BH")
sigPathways$pathway<-rownames(sigPathways)
sigPathways$Program<-paste0(i,"-PosVSNeg")
result<-rbind(result,sigPathways)}

result$pathway=gsub("HALLMARK_","",result$pathway)
result$pathway=str_to_title(result$pathway)



##############
#Plot
#############
k=0
gs=list()
data=data.frame()
for (j in index){
group=paste0(j,"-PosVSNeg")
k=k+1

df1=result[result$Program==group,]
img<-df1[order(df1$t),]
row.names(img)=img$pathway

df2 = data.frame(x = factor(1:nrow(img),labels =row.names(img)), y = img$t,z=img$adj.P.Val)
str(df2)
df2 <- transform(df2, judge = ifelse(y>0, 'Yes', 'No'))
df2$judge=as.character(df2$judge)
df2$judge[df2$z>0.05] = "NS"
df2$judge=factor(df2$judge,levels=c('Yes','No',"NS"))
df2$celltype=j
data=rbind(data,df2)
gs[[k]]=ggplot(data = df2, mapping = aes(x = x, y = y, fill = judge)) + 
  geom_bar(stat = 'identity', position = 'identity') + 
  labs(title=paste0("Epithelial Program: ",group),x="",y="t value of GSVA score")+
  scale_fill_manual(values = c('#fbc430','blue',"gray"), guide = FALSE) +coord_flip()
}
res1=gridExtra::grid.arrange(grobs=gs,ncol=4)
#ggsave('EachEpiProgram_GSVA.pdf',res1,width =8*4 , height = 8*2,limitsize = F)

colnames(data)=c("Hallmark", "t_value", "adj.P.Val", "judge", "celltype")
data$celltype=factor(data$celltype,levels = index)
data1=dcast(data,Hallmark~celltype,value.var = c("t_value"))
rownames(data1)=data1$Hallmark
data1=data1[,-1]
pheatmap(data1,cluster_cols = F,scale = 'column',
            color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                       "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
            fontsize = 12,
            fontsize_col = 12,
            fontsize_row =12,
            cellwidth = 10, 
            cellheight = 10,
            border_color=c("white"),
            treeheight_row=5,
            main="Hallmark pathways"
            
)
#ggsave('Program_GSVA_heatmap.pdf',PP,height = 10,width = 6)
apply(abs(data1),1,range)

#####select pathway####

Epi1=c("Mtorc1_signaling","Oxidative_phosphorylation","Reactive_oxigen_species_pathway","Glycolysis","Pi3k_akt_mtor_signaling")

Epi2=c("Apical_surface","Pi3k_akt_mtor_signaling","Complement","P53_pathway","Apical_junction","Inflammatory_response")

Cycling=c("G2m_checkpoint","E2f_targets","Myc_targets_v1","Mtorc1_signaling","Mitotic_spindle")

Oxd=c("Xenobiotic_metabolism","Reactive_oxigen_species_pathway","Fatty_acid_metabolism","Adipogenesis","Oxidative_phosphorylation")

Stress=c("Tnfa_signaling_via_nfkb","Uv_response_up","P53_pathway","Apoptosis","Uv_response_dn")

AP=c("Allograft_rejection","Interferon_gamma_response","Interferon_alpha_response","Complement","Inflammatory_response")

Innate=c("Apoptosis","Complement","Tnfa_signaling_via_nfkb","Inflammatory_response","P53_pathway")

Mes=c("Epithelial_mesenchymal_transition","Angiogenesis","Coagulation","Wnt_beta_catenin_signaling","Apical_junction")

pathway=unique(c(Stress, Epi1,Oxd, Cycling,Mes, AP, Innate, Epi2))
edit(pathway)
c("Tnfa_signaling_via_nfkb", "Uv_response_up", "P53_pathway", 
  "Apoptosis", "Uv_response_dn", "Mtorc1_signaling", "Oxidative_phosphorylation", 
  "Reactive_oxigen_species_pathway", "Glycolysis", "Pi3k_akt_mtor_signaling", 
  "Xenobiotic_metabolism", "Fatty_acid_metabolism", "Adipogenesis", 
  "G2m_checkpoint", "E2f_targets", "Myc_targets_v1", "Mitotic_spindle", 
  "Epithelial_mesenchymal_transition", "Angiogenesis", "Coagulation", 
  "Wnt_beta_catenin_signaling", "Apical_junction", "Allograft_rejection", 
  "Interferon_gamma_response", "Interferon_alpha_response", "Complement", 
  "Inflammatory_response", "Apical_surface")


pathway1=c("Complement", "Tnfa_signaling_via_nfkb", "Apoptosis","Inflammatory_response", 
            "P53_pathway","Uv_response_up",
           "Allograft_rejection", "Interferon_gamma_response","Interferon_alpha_response",
           "E2f_targets", "G2m_checkpoint", "Myc_targets_v1",
           "Mtorc1_signaling",  "Oxidative_phosphorylation", "Glycolysis", "Pi3k_akt_mtor_signaling",
            "Apical_surface","Apical_junction",
           "Epithelial_mesenchymal_transition",  "Angiogenesis",   "Coagulation", 
           "Xenobiotic_metabolism", "Reactive_oxigen_species_pathway", "Fatty_acid_metabolism", "Adipogenesis")

data2=data1[pathway1,]
pheatmap(data2,scale = 'column',cluster_cols = F,cluster_rows = F,
            color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                       "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
            fontsize = 12,
            fontsize_col = 12,
            fontsize_row =12,
            cellwidth = 10, 
            cellheight = 10,
            border_color=c("white"),
            treeheight_row=5,
            treeheight_col=5,
            main="Hallmark pathways"
            
)
head(data2)
colnames(data2)[colnames(data2)=="Innate"]="Mucosal"
#ggsave('Program_GSVA_heatmap_selected.pdf',P2,height = 10,width = 6)
setwd("/home/data/human/ESCC/Linna-ESCC/Fig/source_data/")
write.table(data2,"Fig3c.txt",quote = F,sep="\t",col.names = NA)

data3=data2
data3[data3>100]=100
pheatmap(data3,
         color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                    "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
         fontsize = 12,
         fontsize_col = 12,
         fontsize_row =12,
         cellwidth = 10, 
         cellheight = 10,
         border_color=c("white"),
         treeheight_row=5,
         treeheight_col=5,
         main="Hallmark pathways"
         
)


