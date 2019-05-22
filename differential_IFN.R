# Looks at differential expression of genes from scRNA-seq data in tabular (gene x cell) format
# Allie Stanton, 2019

#############################################
#       INPUT GENE(S) OF INTEREST HERE      #
#############################################

GOI <- c('IFNAR1')


library(ggplot2)

# read in data
PR.data <- read.table("5835_tagged.dge.txt.gz",header=T,row.names=1)
mock.data <- read.table("5834_tagged.dge.txt.gz",header=T,row.names=1)

# find mitochondrial genes
PR.mito.genes <- grep(pattern = "^MT-", x = rownames(x = PR.data), value = TRUE)
PR.percent.mito <- Matrix::colSums(PR.data[PR.mito.genes, ])/Matrix::colSums(PR.data)

mock.mito.genes <- grep(pattern = "^MT-", x = rownames(x = mock.data), value = TRUE)
mock.percent.mito <- Matrix::colSums(mock.data[mock.mito.genes, ])/Matrix::colSums(mock.data)

# remove cells with more than 2% mitochondrial RNA
PR.data <- PR.data[PR.percent.mito<=0.02]
mock.data <- mock.data[mock.percent.mito<=0.02]

# transpose
PR.tpose <- as.data.frame(t(PR.data))
PR.tpose <- PR.tpose[complete.cases(PR.tpose),]

mock.tpose <- as.data.frame(t(mock.data))
mock.tpose <- mock.tpose[complete.cases(mock.tpose),]

# get only cells scored as positive/negative
PR.pos <- PR.tpose[PR.tpose$GP1>25,]
PR.neg <- PR.tpose[PR.tpose$GP1<=25,]

# create dataframe of cells @ GOI
mock.df <- data.frame("cell"=rownames(mock.tpose), "sample"="mock")
for (gene in GOI){
  mock.df[[gene]] <- mock.tpose[[gene]]
}

PR.pos.df <- data.frame("cell"=rownames(PR.pos), "sample"="PR infected")
for (gene in GOI){
  PR.pos.df[[gene]] <- PR.pos[[gene]]
}

PR.neg.df <- data.frame("cell"=rownames(PR.neg), "sample"="PR uninfected")
for (gene in GOI){
  PR.neg.df[[gene]] <- PR.neg[[gene]]
}

GOI.df <- rbind(mock.df, rbind(PR.pos.df, PR.neg.df))

# box plot based on gene of interest
for (gene in GOI){
  plot <- 
    ggplot(GOI.df, aes(x=GOI.df$sample, y=GOI.df[[gene]])) + 
    geom_boxplot(notch=TRUE) +
    labs(x="sample", y=gene)
      
  print(plot)
}

# bar plot of total copies of gene of interest
for (gene in GOI){
  plot <- 
    ggplot(GOI.df, aes(x=GOI.df$sample, y=GOI.df[[gene]])) + 
    geom_bar(stat="identity") +
    labs(x="sample", y=gene)
  
  print(plot)
}