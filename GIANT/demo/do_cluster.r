# cell classification

# demo START >>>>>>>>>>>>>>>>>>>>>>>>>>
library(Seurat)
library(dplyr)
library(Matrix)

# Load the dataset
#detach("package:gi", unload=TRUE)
library("gi")

projectID <- "geact2"
sessionID <- "5c2325197e71c30a28f5359b"

expr_data <- t(getExprMatrix(projectID, sessionID))
dim(expr_data)
cellStat <- getMeta(projectID, sessionID)
dim(cellStat)

# output dir
#dir.create(path = "03-expression/merged", showWarnings = F, recursive = T)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
expr <- CreateSeuratObject(raw.data = expr_data, min.cells = 3, min.genes = 500, project = "blood", names.delim = "/")
dim(expr@raw.data)

# Standard pre-processing workflow

# QC and selecting cells for further analysis

# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(expr@raw.data))
length(mito.genes)
percent.mito <- Matrix::colSums(expr@raw.data[mito.genes, ])/Matrix::colSums(expr@raw.data)
###

# batch info
batch <- gsub("_cell.*", "", colnames(expr@raw.data))
names(batch) <- colnames(expr@raw.data)
bigBatch <- gsub("-[0-9]*_cell.*", "", colnames(expr@raw.data))
names(bigBatch) <- colnames(expr@raw.data)
#people <- rep(c("A", "B"), c(792, 2463-792))
#names(people) <- colnames(expr@raw.data)
outer_id <- ceiling(as.numeric(gsub(".*cell", "", colnames(expr@raw.data))) / 8); outer_id <- factor(outer_id)
names(outer_id) <- colnames(expr@raw.data)
inner_id <- ceiling(as.numeric(gsub(".*cell", "", colnames(expr@raw.data))) %% 8); inner_id[inner_id==0] <- 8; inner_id <- factor(inner_id, levels = 1:8)
names(inner_id) <- colnames(expr@raw.data)
bigInner_id <- ifelse(as.numeric(inner_id)<=4, "1-4", "5-8")
names(bigInner_id) <- colnames(expr@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
#pdf("03-expression/merged/Seurat_QCstats.pdf",width = 6, height = 6, useDingbats = F)

expr <- AddMetaData(object = expr, metadata = percent.mito, col.name = "percent.mito")
expr <- AddMetaData(object = expr, metadata = batch, col.name = "batch")
expr <- AddMetaData(object = expr, metadata = bigBatch, col.name = "bigBatch")
#expr <- AddMetaData(object = expr, metadata = people, col.name = "people")
expr <- AddMetaData(object = expr, metadata = outer_id, col.name = "outer_id")
expr <- AddMetaData(object = expr, metadata = inner_id, col.name = "inner_id")
expr <- AddMetaData(object = expr, metadata = bigInner_id, col.name = "bigInner_id")
VlnPlot(object = expr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.title.use = 16)

#dev.off()

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
dim(expr@meta.data)
expr <- FilterCells(object = expr, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, Inf))
dim(expr@meta.data)

#Normalizing the data
expr <- NormalizeData(object = expr, normalization.method = "LogNormalize", scale.factor = 10000)

#Detection of variable genes across the single cells
#pdf("03-expression/merged/Seurat_PCA.pdf",width = 6, height = 6, useDingbats = F)

expr <- FindVariableGenes(object = expr, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 1)
length(x = expr@var.genes)

#Scaling the data and removing unwanted sources of variation
#expr <- ScaleData(object = expr, vars.to.regress = c("nUMI", "percent.mito"), num.cores = 20, do.par = T)
expr <- ScaleData(object = expr, vars.to.regress = c("nUMI"), num.cores = 15, do.par = T)

#Perform linear dimensional reduction
expr <- RunPCA(object = expr, pc.genes = expr@var.genes, pcs.compute = 100, do.print = F)

# Examine and visualize PCA results a few different ways
#PrintPCA(object = expr, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
#par(oma=c(0,2,0,0))
#VizPCA(object = expr, pcs.use = 1:2)
#PCAPlot(object = expr, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
expr <- ProjectPCA(object = expr, do.print = FALSE)
#Determine statistically significant principal components
#dev.off()

dims_use <- 1:20

# Cluster the cells
if(! exists("expr_ori")) {
  print("Create copy for original expr")
  expr_ori <- expr
} else {
  print("Use original expr")
  expr <- expr_ori
}
expr <- FindClusters(object = expr, reduction.type = "pca", dims.use = dims_use, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")
#PrintFindClustersParams(object = expr)
expr@meta.data$cluster <- expr@meta.data[, grep("res.", colnames(expr@meta.data), fixed = T)]
#table(expr@meta.data$cluster)

expr <- RunTSNE(object = expr, dims.use = dims_use, do.fast = TRUE)

TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)

# demo END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
