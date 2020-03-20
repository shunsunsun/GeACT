# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("Seurat")
library("dplyr")
library("Matrix")
library("pheatmap")
library("reshape2")
library("grid")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("topGO")
source("../../scripts/cluster_tools.r")
source("../../scripts/pheatmap_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/cellCluster/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/clustering.RData"))

# 1. pre-process ----
# Load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Genomes/human/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(geneID)
colnames(geneID) <- c("ensembl", "symbol")

# Load the dataset
expr_data <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
dim(expr_data)

# Initialize the Seurat object with the raw (non-normalized data).
expr <- CreateSeuratObject(raw.data = expr_data, min.cells = 3, min.genes = 500, project = samplingPos, names.delim = "/")
dim(expr@raw.data)

# Standard pre-processing workflow

# QC and selecting cells for further analysis

# calculate the percent.mito values.
#cellStat <- read.table("03-expression/merged/filtering_cells.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
#dim(cellStat)
#percent.mito <- cellStat[cellStat$filter, "mitoRatio"]
#names(percent.mito) <- rownames(cellStat)[cellStat$filter]

### estimate these two strategies
mito.genes <- grep(pattern = "^MT-", x = rownames(expr@raw.data))
length(mito.genes)
percent.mito <- Matrix::colSums(expr@raw.data[mito.genes, ])/Matrix::colSums(expr@raw.data)
# all(names(percent.mito_alt) == names(percent.mito))
# plot(percent.mito_alt, percent.mito)
# cor(percent.mito_alt, percent.mito)
###

# add metaData
expr <- do_addMeta(expr)

#pdf("03-expression/merged/Seurat_QCstats.pdf",width = 6, height = 6, useDingbats = F)
VlnPlot(object = expr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.title.use = 16)

par(mfrow = c(1, 2))
GenePlot(object = expr, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = expr, gene1 = "nUMI", gene2 = "nGene")
par(mfrow = c(1, 1))

#dev.off()

# filtering cells
dim(expr@meta.data)
expr <- FilterCells(object = expr, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, Inf))
dim(expr@meta.data)

# Normalizing the data
expr <- NormalizeData(object = expr, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
#pdf("03-expression/merged/Seurat_PCA.pdf",width = 6, height = 6, useDingbats = F)

expr <- FindVariableGenes(object = expr, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.25, x.high.cutoff = 5, y.cutoff = 0.5)
length(expr@var.genes)

#Scaling the data and removing unwanted sources of variation
expr <- ScaleData(object = expr, vars.to.regress = c("nUMI", "percent.mito"), num.cores = 20, do.par = T)

#Perform linear dimensional reduction
expr <- RunPCA(object = expr, pc.genes = expr@var.genes, pcs.compute = 100, do.print = F)

# Examine and visualize PCA results a few different ways
PrintPCA(object = expr, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

#par(oma=c(0,2,0,0))
VizPCA(object = expr, pcs.use = 1:2)

PCAPlot(object = expr, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset
expr <- ProjectPCA(object = expr, do.print = FALSE)

# PCHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, 
# and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and genes are ordered according to their PCA scores.
PCHeatmap(object = expr, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = expr, pc.use = 1:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

dmtop <- DimTopGenes(object = expr, dim.use = 1, reduction.type = "pca", num.genes = 30, do.balanced = T)
cat(rev(dmtop[16:30]), sep = "\n")
cat(dmtop[1:15], sep = "\n")

topGene <- sapply(1:50, function(i) {
  y <- DimTopGenes(object = expr, dim.use = i, reduction.type = "pca", num.genes = 30, do.balanced = T)
})
colnames(topGene) <- 1:50
View(topGene)
topGene_DF <- data.frame(PC=rep(1:50, each = 15), dup = duplicated(as.character(topGene)))
topGene_dup <- as.data.frame(table(topGene_DF))
ggplot(topGene_dup, aes(x = PC, y = Freq, fill = dup)) + geom_bar(stat = "identity")

#Determine statistically significant principal components
expr <- JackStraw(object = expr, num.pc = 100, num.replicate = 100, num.cores = 20, do.par = T)
JackStrawPlot(object = expr, PCs = 1:100)
PCElbowPlot(object = expr, num.pc = 100)

cairo_pdf(paste0(OUT, "/Seurat_dim.pdf"), width = 6, height = 10, onefile = T)
JackStrawPlot(object = expr, PCs = 1:100)
PCElbowPlot(object = expr, num.pc = 100) + theme(aspect.ratio = 1)
dev.off()

#dev.off()

save.image(file = paste0(OUT, "/Seurat_step1.RData"))
#load(paste0(OUT, "/Seurat_step1.RData"))

### 2. Cluster the cells ----
if(! exists("expr_ori")) {
  print("Create copy for original expr")
  expr_ori <- expr
}

dims_use <- 1:84  # 1:25
resol <- 0.6
expr <- FindClusters(object = expr_ori, reduction.type = "pca", dims.use = dims_use, 
                     resolution = resol, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")
#PrintFindClustersParams(object = expr)
expr@meta.data$cluster <- expr@meta.data[, grep("res.", colnames(expr@meta.data), fixed = T)]
table(expr@meta.data$cluster)

# Run Non-linear dimensional reduction (tSNE)
expr <- RunTSNE(object = expr, dims.use = dims_use, nthreads = 20, do.fast = T)

TSNEPlot(object = expr, pt.size = 1, do.label = T, no.legend = T, 
         plot.title = paste0("Dimension: 1:", max(dims_use), " Resolution: ", resol))

### manually curation
# xxx <- as.data.frame(expr@dr$tsne@cell.embeddings)
# rm(xxx)
###

# check over-clustering
# library("foreach")
# registerDoSEQ() # avoid connection issue
# expr_xxx <- ValidateClusters(expr, pc.use = dims_use, top.genes = 30, min.connectivity = 0.01, acc.cutoff = 0.85, verbose = TRUE)
# rm(expr_xxx)

# UMAP
#expr <- RunUMAP(object = expr, reduction.use = "pca", dims.use = dims_use, min_dist = 1)
#DimPlot(object = expr, reduction.use = "umap", no.legend = F, do.return = TRUE, 
#        vector.friendly = TRUE, pt.size = 3) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

# add info to meta
cellMetaData <- merge(expr@meta.data, expr@dr$tsne@cell.embeddings, by = 0, sort = F)
cellMetaData$batch <- factor(cellMetaData$batch, levels = unique(cellMetaData$batch))
cellMetaData$cluster <- factor(cellMetaData$cluster, levels = sort(as.numeric(unique(cellMetaData$cluster))))
cellStat <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/filtering_cells.txt"), header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
cellMetaData <- merge(cellMetaData, cellStat[, -c(4, 5)], by.x = 1, by.y = 0, sort = F)
colnames(cellMetaData)[1] <- "cell"
# add tissue/samplingPos/ident
cell_metatable <- read.table("cell_metatable.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cell_metatable)
cellMetaData <- merge(cellMetaData, cell_metatable[, c("cell", "tissue", "samplingPos", "ident")], by = "cell", sort = F)
cellMetaData$tissue <- gsub("_", " ", Hmisc::capitalize(cellMetaData$tissue))
cellMetaData$tissue <- factor(cellMetaData$tissue, levels = sort(unique(cellMetaData$tissue)))

### check batch effect
# cairo_pdf(paste0(OUT, "/Seurat_batchDebug.pdf"), width = 6, height = 6, onefile = T)
# do_batchDebug()
# dev.off()
###

### check marker (expressed ratio)
# marker_genes <- c("PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD14", "CD19")
# do_checkMarker(x = expr_data, marker = marker_genes, cell = "")
###

### find markers for every cluster
# expr.markers_pn <- do_findMarker(expr, only.pos = F, ncpu = 12)
# expr.markers <- subset(expr.markers_pn, avg_logFC > 0)
# table(expr.markers$cluster)
# expr.markers_ftd <- expr.markers[expr.markers$power>=0.4 & expr.markers$avg_logFC>=log(2), ]
# table(expr.markers_ftd$cluster)
# top1 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
# top3 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(3, avg_logFC)
# top5 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(5, avg_logFC)
# top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
###

### cell type specific marker genes
# expr.markers_spec <- do_findSpecMarker(expr, ncpu = 12)
# table(expr.markers_spec$cluster)
# top1_spec <- expr.markers_spec %>% group_by(cluster) %>% top_n(1, avg_logFC)
# top5_spec <- expr.markers_spec %>% group_by(cluster) %>% top_n(5, avg_logFC)
###

### cell-cycle gene list
# data(cc.genes)
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# expr <- CellCycleScoring(object = expr, s.genes = s.genes, g2m.genes = g2m.genes)
###

# the ident of each cluster
cluster_ident_fst5 <- t(sapply(stringr::str_sort(unique(expr@meta.data$cluster), numeric = T), function(x) {
  y0 <- sort(table(subset(cell_metatable, cell %in% rownames(expr@meta.data)[expr@meta.data$cluster == x], "ident", drop = T)), decreasing = T)
  y <- paste0(names(y0), "(", y0, ")")[1:5]
  return(y)
}))
cluster_ident_fst5 <- as.data.frame(cluster_ident_fst5)
View(cluster_ident_fst5)

pdf(paste0(OUT, "/Seurat_tSNE_withCluster.pdf"), width = 6, height = 6, useDingbats = F)

TSNEPlot(object = expr, pt.size = 1, do.label = T, no.legend = T, 
         plot.title = paste0("Dimension: 1:", max(dims_use), " Resolution: ", resol))

dev.off()

ct_color <- c(brewer.pal(name = "Dark2", n = 8), brewer.pal(name = "Paired", n = 9))
names(ct_color) <- levels(cellMetaData$tissue)
ct_color["Kidney"] <- "#F2E100"
ct_color["Large intestine"] <- "tomato"
ct_color["Lung"] <- "#419FDE"
save(ct_color, file = paste0(OUT, "/ct_color.RData"))

pdf(paste0(OUT, "/Seurat_tSNE.pdf"), width = 6, height = 6, useDingbats = F)

# color by tissue
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = tissue)) + geom_point(size = 0.2) +
  theme(legend.position = "bottom", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, legend.box.margin = margin(t = -15), plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color)

dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_simple.pdf"), width = 6, height = 6, useDingbats = F)

# color by tissue
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = tissue)) + geom_point(size = 0.2, show.legend = F) +
  #theme(legend.position = "bottom", legend.justification = c("center")) + 
  theme(aspect.ratio = 1) + 
  #theme(legend.box.margin = margin(t = -15), plot.margin = margin(t = 10, l = -20, r = 20)) + 
  #theme(legend.key.width = unit(0.2, "cm")) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
  #guides(color = guide_legend(ncol = 4, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color) #+ 
  #annotate("segment", x=-Inf,xend=Inf,y=-Inf,yend=-Inf,arrow=arrow(length = unit(0.5, "cm"))) + 
  #annotate("segment", x=-Inf,xend=-Inf,y=-Inf,yend=Inf,arrow=arrow(length = unit(0.5, "cm")))

dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_spMarker.pdf"), width = 6, height = 7, useDingbats = F)

# FeaturePlot
mk_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1", "HBG1")
mk_labels <- c("EPCAM (CD326)", "PECAM1 (CD31)", "PTPRC (CD45)", "COL1A1", "HBG1")
for(i in seq_along(mk_genes)) {
  p <- FeaturePlot_new(object = expr, features.plot = mk_genes[i], new.title = mk_labels[i], cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 1, no.axes = T, do.plot = F, do.return = T)[[1]]
  p <- p + theme(aspect.ratio = 1, plot.title = element_text(size = 40), plot.margin = margin(5,2,-10,2))
  print(p)
}

dev.off()

rm(expr_data, expr_ori)
#save.image(file = paste0(OUT, "/Seurat_step2.RData"))
#load(paste0(OUT, "/Seurat_step2.RData"))

save.image(file = paste0(OUT, "/clustering.RData"))
