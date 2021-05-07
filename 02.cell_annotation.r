library(ggpubr)
library(pheatmap)
source("/home/data/human/Script/Signature/FunctionSignature.R")
source("/home/data/human/ESCC/Linna-ESCC/Fig/Function.R")

pbmc <- readRDS("/home/data/human/ESCC/Linna-ESCC/Fig/Cellinfo_20190924/Cell.rds/Epithelial.rds")
df_cells = FS101_CellsInfos_celltype(pbmc)
TSNEPlot(object = pbmc,do.label = T)

########################
# CALCULATE PCA
########################
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               pcs.compute = 50,
               genes.print = 5)
df_PCA = as.data.frame(pbmc@dr$pca@cell.embeddings)
df_stats = FS102_SampleClusterMatrix(pbmc)
cells = rownames(pbmc@dr$pca@cell.embeddings)
df_PCA$sample = df_cells[cells, "sample"]

#########################################
# DISTANCE BWTWEEN SAMPLES IN PCA SPACE
#########################################

#DISPERSION OF EACH SAMPLE
calculate_Dispersion = function(dataMatrix, cells){
  dataUsed = dataMatrix[,cells]
  center_of_cells = rowMeans(dataUsed)
  df = as.data.frame(apply(dataUsed, 2, function(x){dist(rbind(x, center_of_cells), 
                                                         method = "euclidean")}))
  return(mean(df[,1]))
}

#Bias Of Each Sample...
calculate_Bias = function(dataMatrix, cells){
  dataUsed = dataMatrix[,cells]
  center_of_cells = rowMeans(dataMatrix)
  df = as.data.frame(apply(dataUsed, 2, function(x){dist(rbind(x, center_of_cells), 
                                                         method = "euclidean")}))
  return(mean(df[,1]))
}

df_PCA20 = as.data.frame(t(pbmc@dr$pca@cell.embeddings[,1:30]))

df_PCA_stats = ddply(df_cells, .(sample), summarise, cellcounts = length(cell),
                     bias = calculate_Bias(df_PCA20, cell),
                     dispersion = calculate_Dispersion(df_PCA20, cell))

df_PCA_stats = subset(df_PCA_stats, cellcounts>=100)
df_PCA_stats$dis_to_bias = df_PCA_stats$dispersion/df_PCA_stats$bias
rownames(df_PCA_stats) = df_PCA_stats$sample


Specific_patients = c("P16", "P39", "P130T", "P2", "P127T",
                   "P5","P10","P28","P9","P79", "P76",
                   "P54","P8","P31","P24","P26",
                   "P27","P22","P12","P128T","P21")
Specific_patients = gsub( "T", "",Specific_patients)

df_PCA_stats$TYPE = "Common"
df_PCA_stats[Specific_patients,]$TYPE = "Specific"


######################
# 20210414
######################
df_PCA_stats$intra_cor = df_cor[df_PCA_stats$sample, "intra_cor"]
df_PCA_stats$inter_cor = df_cor[df_PCA_stats$sample, "inter_cor"]
df_PCA_stats$TYPE = factor(df_PCA_stats$TYPE, levels = c("Specific", "Common"))

p1 = ggplot(df_PCA_stats, aes( dispersion,bias, color=TYPE))+
  geom_point(size=2)+
  scale_color_manual(values=c("#facf5a","#4c8492"))

p2 = ggplot(df_PCA_stats, aes(intra_cor,inter_cor, color=TYPE))+
  geom_point(size=2)+
  scale_color_manual(values=c("#facf5a","#4c8492"))

df_PCA_stats$dis_to_bias = df_PCA_stats$bias/df_PCA_stats$dispersion
df_PCA_stats$inter_to_intra = df_PCA_stats$inter_cor/df_PCA_stats$intra_cor


p3 = ggboxplot(df_PCA_stats, 
          x = "TYPE", y = "dis_to_bias",
          fill = "TYPE",
          notch = TRUE)+stat_compare_means()

p4 = ggboxplot(df_PCA_stats, 
          x = "TYPE", y = "inter_to_intra",
          fill = "TYPE",
          notch = TRUE)+stat_compare_means()

grid.arrange(p1, p2, p3, p4, ncol=2)


