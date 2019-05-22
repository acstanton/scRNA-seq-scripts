# read in data
PR.data <- read.table("181210Geh_uganda_dge.dge.txt.gz",header=T,row.names=1)

# set up seurat object
PR <- CreateSeuratObject(raw.data = PR.data, min.cells = 5, min.genes = 200, 
                             project = "181210Geh")

# pre-processing workflow

# find mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = PR@data), value = TRUE)
percent.mito <- Matrix::colSums(PR@raw.data[mito.genes, ])/Matrix::colSums(PR@raw.data)
PR <- AddMetaData(object = PR, metadata = percent.mito, col.name = "percent.mito")

par(mfrow = c(1, 2))
GenePlot(object = PR, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = PR, gene1 = "nUMI", gene2 = "nGene")

# filter to remove cells with too many/few unique genes and too high pct mito genes
PR <- FilterCells(object = PR, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(100, -Inf), high.thresholds = c(5500, 0.02))

# normalize data
PR <- NormalizeData(object = PR, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# find variable genes
PR <- FindVariableGenes(object = PR, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# scale data and regress out variability due to percent mito genes
PR <- ScaleData(object = PR)

# PCA
PR <- RunPCA(object = PR, pc.genes = PR@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
PCAPlot(PR)

# JackStraw
PR <- JackStraw(object = PR, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = PR, PCs = 1:20)

# tSNE
PR <- FindClusters(object=PR,reduction.type="pca",dims.use=1:20,print.output=0,resolution=1,save.SNN=TRUE,force.recalc=T)
PR <- RunTSNE(object=PR,dims.use=1:20,do.fast=T)
TSNEPlot(object=PR)

PR.markers <- FindAllMarkers(object = PR, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
PR.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# violin plot
VlnPlot(object=PR,features.plot="GP1",nCol=6)

# heatmap differential cluster expression
top10 <- PR.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = PR, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# kendall correlation
PR.tpose <- as.data.frame(t(PR.data))
PR.tpose <- PR.tpose[complete.cases(PR.tpose),]

corr <- c()
for (i in 1:ncol(PR.tpose)){
  corr.stats <- cor.test(PR.tpose$GP1,PR.tpose[,i],data=PR.tpose,alternative="two.sided",method="kendall")
  corr <- c(corr,corr.stats$estimate)
}

corr <- as.numeric(corr$x)
corrtable <- data.frame("correlation"=corr)
row.names(corrtable) <- colnames(PR.tpose)
corrtable <- na.omit(corrtable)
corrtable <- corrtable[order(corrtable$correlation), , drop=F]
corrtable <- corrtable[-c(nrow(corrtable)), , drop=F]

# get only ggcells scored as positive/negative

infect.pos <- infect.tpose[infect.tpose$GP1>25,]
infect.neg <- infect.tpose[infect.tpose$GP1<=25,]
