
setwd("/lustre/user/tianf/06-Human_cell_atlas/test_samples")

# cluster ----

library(Seurat)
library(dplyr)
library(Matrix)

# Load the dataset
expr_data <- read.table(file = "03-expression/merged/UMIcount_filtered.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
dim(expr_data)

# Examine the memory savings between regular and sparse matrices
#dense.size <- object.size(x = as.matrix(x = pbmc.data))
#dense.size
#sparse.size <- object.size(x = pbmc.data)
#sparse.size
#dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
expr <- CreateSeuratObject(raw.data = expr_data, min.cells = 3, min.genes = 7000, project = "Test",names.delim = "/")
dim(expr@raw.data)

# Standard pre-processing workflow

# QC and selecting cells for further analysis

# calculate the percent.mito values.
#gene_ID2name <- read.table(file = "Genomes/human/gene_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)
#mito.genes <- gene_ID2name[grep(pattern = "^MT-", x = gene_ID2name$V2), 1]
#mito.genes <- mito.genes[mito.genes%in%rownames(expr@raw.data)]
mito.genes <- grep(pattern = "^MT-", x = rownames(expr@raw.data))
length(mito.genes)
percent.mito <- Matrix::colSums(expr@raw.data[mito.genes, ])/Matrix::colSums(expr@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pdf("03-expression/merged/seurat_QC_stats.pdf",width = 8, height = 6, useDingbats = F)

expr <- AddMetaData(object = expr, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = expr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.title.use = 16)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = expr, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = expr, gene1 = "nUMI", gene2 = "nGene")

dev.off()

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
dim(expr@meta.data)
expr <- FilterCells(object = expr, subset.names = c("nGene", "percent.mito"), low.thresholds = c(7000, -Inf), high.thresholds = c(12000, 0.15))
dim(expr@meta.data)

#Normalizing the data
expr <- NormalizeData(object = expr, normalization.method = "LogNormalize", scale.factor = 10000)

#Detection of variable genes across the single cells
pdf("03-expression/merged/seurat_PCA.pdf",width = 6, height = 6, useDingbats = F)

expr <- FindVariableGenes(object = expr, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 8, y.cutoff = 1.25)
length(x = expr@var.genes)

#Scaling the data and removing unwanted sources of variation
expr <- ScaleData(object = expr, vars.to.regress = c("nUMI", "percent.mito"))

# write meta table
#write.table(x = expr@meta.data, file = "03-expression/merged/seurat_metaData.txt",row.names = T, col.names = NA,quote = F,sep = "\t")

#Perform linear dimensional reduction
expr <- RunPCA(object = expr, pc.genes = expr@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

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
PCHeatmap(object = expr, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = expr, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

#Determine statistically significant principal components

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
expr <- JackStraw(object = expr, num.replicate = 100, do.print = T)

# The JackStrawPlot function provides a visualization tool
# for comparing the distribution of p-values for each PC with a uniform distribution (dashed line).
# ‘Significant’ PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line).
JackStrawPlot(object = expr, PCs = 1:15)

# A more ad hoc method for determining which PCs to use is to look at a plot
# of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph.
PCElbowPlot(object = expr)

dev.off()

# Cluster the cells
# The FindClusters function implements the procedure, and contains a resolution parameter 
# that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters.
# We find that setting this parameter between 0.6-1.2 typically returns good results
# for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
# The clusters are saved in the object@ident slot.

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
expr <- FindClusters(object = expr, reduction.type = "pca", dims.use = 1:11, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")

PrintFindClustersParams(object = expr)

# Run Non-linear dimensional reduction (tSNE)
# Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets.
# While we no longer advise clustering directly on tSNE components,
# cells within the graph-based clusters determined above should co-localize on the tSNE plot.
# This is because the tSNE aims to place cells with similar local neighborhoods
# in high-dimensional space together in low-dimensional space. As input to the tSNE,
# we suggest using the same PCs as input to the clustering analysis,
# although computing the tSNE based on scaled gene expression is also supported using the genes.use argument.

expr <- RunTSNE(object = expr, dims.use = 1:11, do.fast = TRUE)

pdf("03-expression/merged/seurat_tSNE.pdf",width = 6, height = 6, useDingbats = F)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = expr, pt.size = 2)

# You can save the object at this point.
# save(expr, file = "03-expression/merged/seurat_expr.Robj")

# Finding differentially expressed genes (cluster biomarkers)

# find all markers of cluster 1
#cluster1.markers <- FindMarkers(object = expr, ident.1 = 1, min.pct = 0.25)
#print(x = head(x = cluster1.markers, n = 5))
# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(object = expr, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
expr.markers <- FindAllMarkers(object = expr, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
dim(expr.markers)
expr.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
# filter by q-value
expr.markers_ftd <- expr.markers[expr.markers$p_val_adj<=0.05,]
dim(expr.markers_ftd)

top1 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
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
top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = expr, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# known cell-cycle gene list from literature
s.genes <- read.table("Data/human/cell_cycle_gene_G1_S.txt", header = F, sep = "\t", stringsAsFactors = F)[,1]
g2m.genes <- read.table("Data/human/cell_cycle_gene_G2_M.txt", header = F, sep = "\t", stringsAsFactors = F)[,1]

DoHeatmap(object = expr, genes.use = c(s.genes, g2m.genes), slim.col.label = TRUE, remove.key = TRUE, cex.row = 4)

# Assigning cell type identity to clusters
current.cluster.ids <- c(0, 1)
new.cluster.ids <- c("G2/M", "G1/S")
expr@ident <- plyr::mapvalues(x = expr@ident, from = current.cluster.ids, to = new.cluster.ids)
expr@ident <- factor(expr@ident, levels = c("G1/S", "G2/M"))
expr@meta.data$ident <- expr@ident
TSNEPlot(object = expr, pt.size = 2)
expr.markers_ftd$cluster <- plyr::mapvalues(x = expr.markers_ftd$cluster, from = current.cluster.ids, to = new.cluster.ids)
expr.markers_ftd$cluster <- factor(expr.markers_ftd$cluster, levels = c("G1/S", "G2/M"))
expr.markers_ftd <- expr.markers_ftd[order(expr.markers_ftd$cluster), ]
DoHeatmap(object = expr, genes.use = c(s.genes, g2m.genes), slim.col.label = TRUE, remove.key = TRUE, cex.row = 4)

# compare markers with literature
expr.markers_plus <- unique(c(expr.markers_ftd$gene,s.genes,g2m.genes))
expr.markers_in_literature <- apply(as.matrix(expr.markers_plus), 1, function(x) {
  if(x[1]%in%s.genes) { y <- "G1/S" } else if(x[1]%in%g2m.genes) { y <- "G2/M" } else { y <- "NA" }
  return(y)
})
expr.markers_in_work <- apply(as.matrix(expr.markers_plus), 1, function(x) {
  if(! x[1]%in%expr.markers_ftd$gene) { y <- "NA" } else { y <- as.character(expr.markers_ftd$cluster)[expr.markers_ftd$gene==x[1]] }
  return(y)
})
expr.markers_DF <- data.frame(row.names = expr.markers_plus, literature=expr.markers_in_literature, work=expr.markers_in_work, stringsAsFactors = F)
expr.markers_DF$literature <- factor(expr.markers_DF$literature, levels = c("G1/S","G2/M","NA"))
expr.markers_DF$work <- factor(expr.markers_DF$work, levels = c("G1/S","G2/M","NA"))
expr.markers_TB <- table(expr.markers_DF)
expr.markers_TB

# in literature
top1_in_literature <- expr.markers_ftd[expr.markers_ftd$gene%in%c(s.genes, g2m.genes), ] %>% group_by(cluster) %>% top_n(1, avg_logFC)
VlnPlot(object = expr, features.plot = top1_in_literature$gene, size.title.use = 16)
for(i in top1_in_literature$gene) {
  FeaturePlot(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2)
}

# not in literature
top1_notIn_literature <- expr.markers_ftd[ ! expr.markers_ftd$gene%in%c(s.genes, g2m.genes), ] %>% group_by(cluster) %>% top_n(1, avg_logFC)
VlnPlot(object = expr, features.plot = top1_notIn_literature$gene, size.title.use = 16)
for(i in top1_notIn_literature$gene) {
  FeaturePlot(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2)
}

dev.off()

# save and write meta table
save(expr, file = "03-expression/merged/seurat_expr.Robj")
write.table(x = expr@meta.data, file = "03-expression/merged/seurat_metaData.txt",row.names = T, col.names = NA,quote = F,sep = "\t")
write.table(x = expr.markers_ftd, file = "03-expression/merged/marker_genes.txt",  row.names = F, col.names = T, quote = F, sep = "\t")









# Assigning cell type identity to clusters

###
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)
###

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
