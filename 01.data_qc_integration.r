args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(dplyr)
samplename = args[1]

sourceDir = paste("/home/data/human/data/", samplename,"/GRCh38", sep="")
resultDir = paste("/home/data/human/seurat/", samplename, sep="")

#Check if the folder exists
if (file.exists(resultDir)){
  setwd(resultDir)
} else {
  dir.create(resultDir, showWarnings = FALSE)
  setwd(resultDir)
}

#Check if the RDS exists

pbmc.data <- Read10X(data.dir = sourceDir)

dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size

sparse.size <- object.size(x = pbmc.data)
sparse.size

dense.size/sparse.size

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 500, 
                           project = "10X_PBMC")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

pdf("QC_vln.pdf")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

pdf("QC_MITO.pdf")
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
dev.off()

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.5))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pdf("QC_VARgenes.pdf")
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)    
dev.off()

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

pdf("PCA.pdf")
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
dev.off()


pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = pbmc, PCs = 1:12)

pdf("PCElbowPlot.pdf")
PCElbowPlot(object = pbmc)
dev.off()

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

pdf("TSNEPlot.pdf")
TSNEPlot(object = pbmc,pt.size = 0.5)
dev.off()

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, 
                               min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

png("HEATMAP.png", width = 1000, height = 1200)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

saveRDS(pbmc, file = "Sample.rds")
