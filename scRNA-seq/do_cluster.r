# cell classification
setwd("~/lustre/06-Human_cell_atlas/test_FACS/")

# cluster ----
library(Seurat)
library(dplyr)
library(Matrix)

# Load the dataset
expr_data <- read.table(file = "03-expression/merged/UMIcount_filtered.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
dim(expr_data)
#colnames(expr_data) <- gsub(pattern = "cell", replacement = "_cell", x = colnames(expr_data))

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
expr <- CreateSeuratObject(raw.data = expr_data, min.cells = 3, min.genes = 500, project = "test_FACS", names.delim = "/")
dim(expr@raw.data)

# Standard pre-processing workflow

# QC and selecting cells for further analysis

# calculate the percent.mito values.
cellStat <- read.table("03-expression/merged/filtering_cells.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
percent.mito <- cellStat[cellStat$filter, "mito"]
names(percent.mito) <- rownames(cellStat)[cellStat$filter]

### estimate these two strategies
# mito.genes <- grep(pattern = "^MT-", x = rownames(expr@raw.data))
# length(mito.genes)
# percent.mito_alt <- Matrix::colSums(expr@raw.data[mito.genes, ])/Matrix::colSums(expr@raw.data)
# all(names(percent.mito_alt) == names(percent.mito))
# plot(percent.mito_alt, percent.mito)
# cor(percent.mito_alt, percent.mito)
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

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = expr, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = expr, gene1 = "nUMI", gene2 = "nGene")
par(mfrow = c(1, 1))

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
                          x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 1)
length(x = expr@var.genes)

#Scaling the data and removing unwanted sources of variation
expr <- ScaleData(object = expr, vars.to.regress = c("nUMI", "percent.mito"), num.cores = 20, do.par = T)

#Perform linear dimensional reduction
expr <- RunPCA(object = expr, pc.genes = expr@var.genes, pcs.compute = 100, do.print = F)

# Examine and visualize PCA results a few different ways
PrintPCA(object = expr, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

#par(oma=c(0,2,0,0))
VizPCA(object = expr, pcs.use = 1:2)

PCAPlot(object = expr, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
expr <- ProjectPCA(object = expr, do.print = FALSE)

# PCHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, 
# and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and genes are ordered according to their PCA scores.
PCHeatmap(object = expr, pc.use = 1, cells.use = 232, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = expr, pc.use = 1:12, cells.use = 232, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

#Determine statistically significant principal components

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
expr <- JackStraw(object = expr, num.replicate = 100, num.cores = 20, do.par = T, maxit = 1000)

# The JackStrawPlot function provides a visualization tool
# for comparing the distribution of p-values for each PC with a uniform distribution (dashed line).
# ‘Significant’ PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line).
JackStrawPlot(object = expr, PCs = 1:20)

# A more ad hoc method for determining which PCs to use is to look at a plot
# of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph.
PCElbowPlot(object = expr, num.pc = 26)

#dev.off()

dims_use <- 1:25

# Cluster the cells
# The FindClusters function implements the procedure, and contains a resolution parameter 
# that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters.
# We find that setting this parameter between 0.6-1.2 typically returns good results
# for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
# The clusters are saved in the object@ident slot.

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
if(! exists("expr_ori")) {
  print("Create copy for original expr")
  expr_ori <- expr
} else {
  print("Use original expr")
  expr <- expr_ori
}
expr <- FindClusters(object = expr, reduction.type = "pca", dims.use = dims_use, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")
PrintFindClustersParams(object = expr)
table(expr@meta.data$res.0.6)

# Run Non-linear dimensional reduction (tSNE)
# Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets.
# While we no longer advise clustering directly on tSNE components,
# cells within the graph-based clusters determined above should co-localize on the tSNE plot.
# This is because the tSNE aims to place cells with similar local neighborhoods
# in high-dimensional space together in low-dimensional space. As input to the tSNE,
# we suggest using the same PCs as input to the clustering analysis,
# although computing the tSNE based on scaled gene expression is also supported using the genes.use argument.

expr <- RunTSNE(object = expr, dims.use = dims_use, do.fast = TRUE)
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)

### check batch effect
cellMetaData <- merge(expr@meta.data, expr@dr$tsne@cell.embeddings, by = 0, sort = F)
#cellMetaData$batch <- factor(cellMetaData$batch, levels = unique(cellMetaData$batch))
cellStat <- read.table("03-expression/merged/filtering_cells.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
cellMetaData <- merge(cellMetaData, cellStat[, -c(4, 5)], by.x = 1, by.y = 0, sort = F)

pdf("03-expression/merged/Seurat_batchEffect.pdf", width = 6, height = 6, useDingbats = F)
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = res.0.6)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = batch)) + geom_point(show.legend = T) + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect(fill=alpha('white', 0.4))) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = bigBatch)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
#ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = people)) + geom_point() + theme_bw() + 
#  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
#  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = outer_id)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect(fill = alpha('white', 0.4))) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = inner_id)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect(fill = alpha('white', 0.4))) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = bigInner_id)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect(fill = NA)) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))

# expression info
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = nGene)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top"))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = nUMI)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top"))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = ABratio)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top"))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = cleanReads)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top"))
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = mpRatio)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect()) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top"))

dev.off()
###

### check marker (expressed ratio)
do_checkMarker <- function(x, marker, cell) {
  expr_sub <- x[marker, grep(cell, colnames(x))]
  rownames(expr_sub) <- marker
  #cat("Data dimension: ", dim(expr_sub), "\n")
  expr_ratio <- t(apply(expr_sub, 1, function(x) { y0 <- sum(x==0); y1 <- sum(x>0); y2 <- y1/(y0 + y1); return(c(y0, y1, y2)) }))
  colnames(expr_ratio) <- c("Non_expressed", "Expressed", "ratio")
  expr_ratio <- data.frame(bigBatch=cell, marker=rownames(expr_ratio), expr_ratio)
  return(expr_ratio)
}

marker_genes <- c("PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD14", "CD19")
do_checkMarker(x = expr_data, marker = marker_genes, cell = "FACS")
###

pdf("03-expression/merged/Seurat_tSNE.pdf",width = 6, height = 6, useDingbats = F)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)
table(expr@meta.data$res.0.6)

# Finding differentially expressed genes (cluster biomarkers)
# find all markers of cluster 1
#cluster1.markers <- FindMarkers(object = expr, ident.1 = 1, min.pct = 0.25)
#print(x = head(x = cluster1.markers, n = 5))
# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(object = expr, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report only the positive ones
expr.markers <- FindAllMarkers(object = expr, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
dim(expr.markers)
table(expr.markers$cluster)
expr.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
# filter by q-value
#expr.markers_ftd <- expr.markers[expr.markers$p_val_adj<=0.05,]
#dim(expr.markers_ftd)
#table(expr.markers_ftd$cluster)
# filter by power
expr.markers_ftd <- expr.markers[expr.markers$power>=0.4 & expr.markers$avg_logFC>=log(1.5), ]
dim(expr.markers_ftd)
table(expr.markers_ftd$cluster)

top1 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
# Seurat has four tests for differential expression which can be set with the test.use parameter:
# ROC test (“roc”), t-test (“t”), LRT test based on zero-inflated data (“bimod”, default),
# LRT test based on tobit-censoring models (“tobit”) The ROC test returns the ‘classification power’
# for any individual marker (ranging from 0 - random, to 1 - perfect).
#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)

# We include several tools for visualizing marker expression.
# VlnPlot (shows expression probability distributions across clusters),
# and FeaturePlot (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations.
# We also suggest exploring JoyPlot, CellPlot, and DotPlot as additional methods to view your dataset.
VlnPlot(object = expr, features.plot = top1$gene, size.title.use = 16)
# you can plot raw UMI counts as well
VlnPlot(object = expr, features.plot = top1$gene, use.raw = TRUE, y.log = TRUE, size.title.use = 16)

for(i in top1$gene) {
  FeaturePlot(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2)
}

# DoHeatmap generates an expression heatmap for given cells and genes.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
#top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = expr, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = expr, genes.use = c(top10$gene, "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA"), slim.col.label = TRUE, remove.key = T)
#DoHeatmap(object = expr, genes.use = c(expr.markers_ftd$gene, "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA"), slim.col.label = TRUE, remove.key = T)

## find marker genes for sub-types
# markers_1_vs_2_3 <- FindMarkers(object = expr, test.use = "roc", only.pos = TRUE, ident.1 = 1, ident.2 = c(2, 3), min.pct = 0.25)
# subset(markers_1_vs_2_3, power>=0.4 & avg_logFC>=log(2))
# markers_2_vs_1_3 <- FindMarkers(object = expr, test.use = "roc", only.pos = TRUE, ident.1 = 2, ident.2 = c(1, 3), min.pct = 0.25)
# subset(markers_2_vs_1_3, power>=0.4 & avg_logFC>=log(2))
# markers_3_vs_1_2 <- FindMarkers(object = expr, test.use = "roc", only.pos = TRUE, ident.1 = 3, ident.2 = c(1, 2), min.pct = 0.25)
# subset(markers_3_vs_1_2, power>=0.4 & avg_logFC>=log(2))

### known cell-cycle gene list from literature
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DoHeatmap(object = expr, genes.use = c(s.genes, g2m.genes), slim.col.label = TRUE, remove.key = TRUE, cex.row = 4)

dev.off()

### check abnormal sub-types
pdf("03-expression/merged/Seurat_mito.pdf",width = 6, height = 6, useDingbats = F)

ggplot(cellMetaData, aes(x = res.0.6, y = percent.mito, fill = res.0.6)) + geom_boxplot(show.legend = F, notch = T)
ggplot(cellMetaData, aes(x = res.0.6, y = nUMI, fill = res.0.6)) + geom_boxplot(show.legend = F, notch = T)
ggplot(cellMetaData, aes(x = res.0.6, y = nUMI*percent.mito, fill = res.0.6)) + geom_boxplot(show.legend = F, notch = T)

# library(gtools)
# tmp <- table(gsub("cell.*", "", WhichCells(expr, ident = 0)))
# plot(tmp[mixedsort(names(tmp))], las = 2, ylab = "Cell number")

dev.off()


pdf(file = "03-expression/merged/Seurat_clusterDebug.pdf", width = 5, height = 4, useDingbats = F)
## bigBatch level
# outer
for(i in unique(cellMetaData$bigBatch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
    geom_violin(aes(x = outer_id, y = mito), fill = "navy", alpha = 0.3) + 
    ggtitle(paste(i, "(outer barcode: mito)"))
  print(gp)
}

for(i in unique(cellMetaData$bigBatch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
    geom_bar(aes(x = outer_id, fill = res.0.6)) + 
    ggtitle(paste(i, "(outer barcode)")) + 
    scale_y_continuous(expand = c(0, 0)) + 
    ylab("Cell number")
  print(gp)
}

for(i in unique(cellMetaData$bigBatch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
    geom_bar(aes(x = outer_id, fill = res.0.6), position="fill") + 
    ggtitle(paste(i, "(outer barcode)")) + 
    scale_y_continuous(expand = c(0, 0)) + 
    ylab("Cell frequency")
  print(gp)
}

# inner
for(i in unique(cellMetaData$bigBatch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
    geom_violin(aes(x = inner_id, y = mito), fill = "navy", alpha = 0.3) + 
    ggtitle(paste(i, "(inner barcode: mito)"))
  print(gp)
}

for(i in unique(cellMetaData$bigBatch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
    geom_bar(aes(x = inner_id, fill = res.0.6)) + 
    ggtitle(paste(i, "(inner barcode)")) + 
    scale_y_continuous(expand = c(0, 0)) + 
    ylab("Cell number")
  print(gp)
}

for(i in unique(cellMetaData$bigBatch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
    geom_bar(aes(x = inner_id, fill = res.0.6), position="fill") + 
    ggtitle(paste(i, "(inner barcode)")) + 
    scale_y_continuous(expand = c(0, 0)) + 
    ylab("Cell frequency")
  print(gp)
}

## batch level
for(i in unique(cellMetaData$batch)) {
  print(i)
  gp <- ggplot(cellMetaData[cellMetaData$batch==i, ]) + 
    geom_bar(aes(x = 1, y = sqrt(cleanReads), fill = mito, color = res.0.6), alpha = 1, stat = "identity", show.legend = T) + 
    coord_polar("x") + 
    geom_text(x=0, y=0, vjust = 2, aes(label = round(mito*100,2)), size=2) + 
    facet_grid(inner_id ~ outer_id, drop = F) + ggtitle(i) + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
    theme(strip.text = element_text(size = 6)) + 
    theme(plot.margin = unit(c(0.1,0.1,0,0), "cm")) + 
    scale_fill_gradientn(colours = c("white", "blue"))
  print(gp)
}

dev.off()

### Assigning cell type identity to clusters
if(! exists("expr_unassigned")) {
  print("Create copy for unassigned expr")
  expr_unassigned <- expr
}

# check marker genes
expr.markers_ftd %>% group_by(cluster) %>% top_n(2, avg_logFC)

library("enrichR")
#dbs <- listEnrichrDbs()
enriched_LS <- lapply(split(x = expr.markers_ftd$gene, f = expr.markers_ftd$cluster), function(x) {
  y <- enrichr(genes = x, databases = "GO_Biological_Process_2018")
  return(y)
})
lapply(enriched_LS, function(x) { x[[1]][1:15, 1:4] })
###

### marker genes in similar clusters
# # 1
# tmp <- FindMarkers(object = expr_unassigned, test.use = "roc", only.pos = TRUE, min.pct = 0.25, ident.1 = 1, ident.2 = 0)
# tmp_ftd <- subset(tmp, power>=0.4 & avg_logFC>=log(2))
# dim(tmp_ftd)
# FeaturePlot(object = expr_unassigned, features.plot = "UTP18", cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, 
#             cells.use = WhichCells(object = expr_unassigned, ident = c(0, 3)))
# 
# tmp <- FindMarkers(object = expr_unassigned, test.use = "roc", only.pos = TRUE, min.pct = 0.25, ident.1 = 3, ident.2 = 0)
# tmp_ftd <- subset(tmp, power>=0.4 & avg_logFC>=log(2))
# dim(tmp_ftd)
# FeaturePlot(object = expr_unassigned, features.plot = "ANXA1", cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, 
#             cells.use = WhichCells(object = expr_unassigned))
# # [merge]
# tmp <- FindMarkers(object = expr_unassigned, test.use = "roc", only.pos = TRUE, min.pct = 0.25, ident.1 = c(1,4), ident.2 = 6)
# tmp_ftd <- subset(tmp, power>=0.4 & avg_logFC>=log(2))
# FeaturePlot(object = expr_unassigned, features.plot = "CCR7", cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, 
#             cells.use = WhichCells(object = expr_unassigned, ident = c(1,4,6)))
###

#levels(expr@ident) <- seq(0, length(levels(expr@ident)) - 1)
#expr <- expr_unassigned
current.cluster.ids <- as.numeric(levels(expr@ident))
new.cluster.ids <- c("Test")
expr@ident <- plyr::mapvalues(x = expr@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = expr, do.label = TRUE, pt.size = 1, no.legend = T)

cellMetaData <- merge(cellMetaData, data.frame(expr@ident), by.x = 1, by.y = 0, sort = F)
expr.markers_ftd_labelled <- merge(x = expr.markers_ftd, y = unique(cellMetaData[, c("res.0.6", "expr.ident")]), by.x = "cluster", by.y = "res.0.6", sort = F)

# save and write meta table
save(expr, file = "03-expression/merged/Seurat_expr.Robj")
write.table(x = cellMetaData, file = "03-expression/merged/Seurat_metaData.txt",row.names = F, col.names = T, quote = F,sep = "\t")
write.table(x = expr.markers_ftd_labelled, file = "03-expression/merged/Seurat_markerGenes.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = expr@scale.data, file = "03-expression/merged/UMIcount_scaled.txt", row.names = T, col.names = NA, quote = F, sep = "\t")

pdf("03-expression/merged/Seurat_tSNE_assigned.pdf",width = 6, height = 6, useDingbats = F)
TSNEPlot(object = expr, do.label = TRUE, pt.size = 1, no.legend = T)
dev.off()

################################################################

####
do_future <- function() {

# Further subdivisions within cell types
# If you perturb some of our parameter choices above
# (for example, setting resolution=0.8 or changing the number of PCs),
# you might see the CD4 T cells subdivide into two groups.
# You can explore this subdivision to find markers separating the two T cell subsets.
# However, before reclustering (which will overwrite object@ident),
# we can stash our renamed identities to be easily recovered later.

# First lets stash our identities for later
pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.8, print.output = FALSE)

# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", 
                  no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)

# Find discriminating markers
tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
FeaturePlot(object = pbmc, features.plot = c("S100A4", "CCR7"), cols.use = c("green", "blue"))

# The memory/naive split is bit weak, and we would probably benefit from looking at more cells
# to see if this becomes more convincing. In the meantime,
# we can restore our old cluster identities for downstream processing.
pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
save(pbmc, file = "Seurat_tutorial/pbmc3k_final.Rda")

# Visualization ----
# http://satijalab.org/seurat/visualization_vignette.html

pbmc

features.plot <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
# Joy plots - from ggjoy. Visualize single cell expression distributions in
# each cluster
JoyPlot(object = pbmc, features.plot = features.plot, nCol = 2)

# Violin plots. Visualize single cell expression distributions in each
# cluster
VlnPlot(object = pbmc, features.plot = features.plot, x.lab.rot = TRUE)

# Dot plots - the size of the dot corresponds to the percentage of cells
# expressing the gene in each cluster. The color represents the average
# expression level
DotPlot(object = pbmc, genes.plot = features.plot, plot.legend = TRUE)

# Feature plot - visualize gene expression in low-dimensional space
FeaturePlot(object = pbmc, features.plot = features.plot, cols.use = c("lightgrey", "blue"))

# Single cell heatmap of gene expression
DoHeatmap(object = SubsetData(object = pbmc, max.cells.per.ident = 100), genes.use = features.plot, slim.col.label = TRUE, group.label.rot = TRUE)

## New additions to FeaturePlot

# Plot a legend to map colors to expression levels
FeaturePlot(object = pbmc, features.plot = "MS4A1", no.legend = FALSE)

# Adjust the contrast in the plot
FeaturePlot(object = pbmc, features.plot = "MS4A1", no.legend = FALSE, min.cutoff = 1, 
            max.cutoff = 3)

# Calculate gene-specific contrast levels based on quantiles of non-zero
# expression. Particularly useful when plotting multiple markers
FeaturePlot(object = pbmc, features.plot = c("MS4A1", "PTPRCAP"), no.legend = FALSE, 
            min.cutoff = "q10", max.cutoff = "q90")

# Visualize co-expression of two genes simultaneously. in beta mode (data is
# currently discretized, will be placed on a continuous scale soon)
FeaturePlot(object = pbmc, features.plot = c("MS4A1", "CD79A"), cols.use = c("grey", 
                                                                             "red", "blue", "green"), overlay = TRUE, no.legend = FALSE)
# Inspire fear and awe at your next lab meeting.
FeaturePlot(object = pbmc, features.plot = c("MS4A1", "PTPRCAP"), no.legend = FALSE, 
            min.cutoff = "q10", max.cutoff = "q90", dark.theme = TRUE)
# Dark themes are also available with the `DarkTheme` function, and can be
# added to any ggplot2-based plot.

## Interactive plotting features
# Seurat utilizes R’s plotly graphing library to create interactive plots.
# Just pass do.hover = TRUE to FeaturePlot, CellPlot,  GenePlot, or DimPlot
# and its extensions (includes PCAPlot, TSNEPlot, etc.)

# Include additional data to display alongside cell names in data.hover
FeaturePlot(object = pbmc, features.plot = "MS4A1", do.hover = TRUE, data.hover = c("ident", "PC1", "nGene"))

# Manually select cells for further investigation.
# We have found this particularly useful for small clusters
# that do not always separate using unbiased clustering, but which look tantalizingly distinct.
# You can now select these cells by passing  do.identify = TRUE to FeaturePlot
# or DimPlot and its extensions (includes PCAPlot, TSNEPlot, etc.).
# You can then set them to a new identity class and perform differential expression.
pbmc <- RenameIdent(object = pbmc, old.ident.name = "Dendritic cells", new.ident.name = "CD14+ Monocytes")
select.cells <- TSNEPlot(object = pbmc, do.identify = TRUE)

###
# manually select
###

head(select.cells)
pbmc <- SetIdent(object = pbmc, cells.use = select.cells, ident.use = "NewCells")

# Now, we find markers that are specific to the new cells, and find clear DC
# markers
newcells.markers <- FindMarkers(object = pbmc, ident.1 = "NewCells", ident.2 = "CD14+ Monocytes", 
                                min.diff.pct = 0.3, only.pos = TRUE)
head(x = newcells.markers)

TSNEPlot(object = pbmc, do.return = T, no.legend = F, do.label = TRUE)

# Dimensional Reduction ----
# http://satijalab.org/seurat/dim_reduction_vignette.html
head(x = GetCellEmbeddings(object = pbmc, reduction.type = "pca", dims.use = 1:5))

head(x = GetGeneLoadings(object = pbmc, reduction.type = "pca", dims.use = 1:5))

# We also provide shortcut functions for common dimensional reduction
# techniques like PCA PCAEmbed and PCALoad() will pull the PCA cell
# embeddings and gene loadings respectively
head(x = GetDimReduction(object = pbmc, reduction.type = "pca", slot = "sdev"))

# Seurat provides RunPCA (pca), RunICA (ica), RunTSNE (tsne), and RunDiffusionMap (dmap),
# representing dimensional reduction techniques commonly applied to scRNA-seq data.
# When using these functions, all slots are filled automatically.

# We also allow users to add the results of a custom dimensional reduction technique
# (for example, multi-dimensional scaling (MDS), or zero-inflated factor analysis),
# that is computed separately. All you need is a matrix with each cell’s coordinates
# in low-dimensional space, as shown below.

## Storing a new dimensional reduction calculation

# MDS
# Before running MDS, we first calculate a distance matrix between all pairs
# of cells.  Here we use a simple euclidean distance metric on all genes,
# using object@scale.data as input
d <- dist(x = t(x = pbmc@scale.data))
# Run the MDS procedure, k determines the number of dimensions
mds <- cmdscale(d = d, k = 2)
# cmdscale returns the cell embeddings, we first label the columns to ensure
# downstream consistency
colnames(x = mds) <- paste0("MDS", 1:2)
# We will now store this as a new dimensional reduction called 'mds'
pbmc <- SetDimReduction(object = pbmc, reduction.type = "mds", slot = "cell.embeddings", 
                        new.data = mds)
pbmc <- SetDimReduction(object = pbmc, reduction.type = "mds", slot = "key", 
                        new.data = "MDS")

# We can now use this as you would any other dimensional reduction in all
# downstream functions (similar to PCAPlot, but generalized for any
# reduction)
DimPlot(object = pbmc, reduction.use = "mds", pt.size = 0.5)

# If you wold like to observe genes that are strongly correlated with the
# first MDS coordinate (similar to ProjectPCA, but generalized for any
# reduction):
pbmc <- ProjectDim(object = pbmc, reduction.type = "mds")

# Display the results as a heatmap (similar to PCHeatmap, but generalized
# for any dimensional reduction)
DimHeatmap(object = pbmc, reduction.type = "mds", dim.use = 1, cells.use = 500, 
           use.full = TRUE, do.balanced = TRUE, label.columns = FALSE, remove.key = TRUE)

# Explore how the first MDS dimension is distributed across clusters
VlnPlot(object = pbmc, features.plot = "MDS1", x.lab.rot = TRUE)

# See how the first MDS dimension is correlated with the first PC dimension
GenePlot(object = pbmc, gene1 = "MDS1", gene2 = "PC1")

## Changes to PCA
# For large datasets containing rare cell types, we often see improved results
# by setting this to FALSE, as this prevents the initial PCs
# (which often explain a disproportionate amount of variance)
# from masking rare cell types or subtle sources of heterogeneity that appear in later PCs.

}
####
