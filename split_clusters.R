# Clusters cells from scRNA-seq data and then splits the original data up into different tables based on the clusters
# Allie Stanton, 2019

# read in data
uganda.data <- read.table("181210Geh_uganda_dge.dge.txt.gz",header=T,row.names=1)

# set up seurat object
uganda <- CreateSeuratObject(counts = uganda.data, min.cells = 5, min.features = 200)

# pre-processing workflow

# find mitochondrial genes
uganda[["percent.mt"]] <- PercentageFeatureSet(uganda, pattern = "^MT-")

# filter to remove cells with too many/few unique genes and too high pct mito genes
uganda <- subset(uganda, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 2)

# normalize data
uganda <- NormalizeData(object = uganda, normalization.method = "LogNormalize", 
                    scale.factor = 10000)

# find variable genes
uganda <- FindVariableFeatures(uganda, selection.method = "vst", nfeatures = 2000)

# scale data and regress out variability due to percent mito genes
all.genes <- rownames(uganda)
uganda <- ScaleData(uganda, features = all.genes)

# PCA
uganda <- RunPCA(uganda, features = VariableFeatures(object = uganda))

# clustering
uganda <- FindNeighbors(uganda, dims = 1:10)
uganda <- FindClusters(uganda, resolution = 1)

uganda <- RunTSNE(uganda, dims = 1:10)
TSNEPlot(uganda)

# find cells in each cluster
cells.by.cluster <- as.data.frame(Idents(uganda))
cells.by.cluster <- split(cells.by.cluster, cells.by.cluster$`Idents(uganda)`)

labels <- paste("cluster", seq_along(cells.by.cluster), sep=".")
for (i in seq_along(cells.by.cluster)){
  cluster.data <- uganda.data[,rownames(as.data.frame(cells.by.cluster[i]))]
  assign(labels[i], cluster.data)
}

# output will be one table for each cluster that you can then save, turn into separate Seurat objects, etc.

