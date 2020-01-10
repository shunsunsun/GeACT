# cell classification
setwd("~/lustre/06-Human_cell_atlas/datasets/09_D1_45/")

library("Seurat")
library("dplyr")
library("Matrix")
library("pheatmap")
source("../../scripts/cluster_tools.r")

samplingPos <- "SD"
OUT <- paste0("03-expression/merged/cellCluster/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

# 1. pre-process ----
# Load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Data/gene_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)
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
expr <- JackStraw(object = expr, num.pc = 30, num.replicate = 100, num.cores = 20, do.par = T)
JackStrawPlot(object = expr, PCs = 1:30)
PCElbowPlot(object = expr, num.pc = 30)

#dev.off()

save.image(file = paste0(OUT, "/Seurat_step1.RData"))
#load(paste0(OUT, "/Seurat_step1.RData"))

### 2. Cluster the cells ----
if(! exists("expr_ori")) {
  print("Create copy for original expr")
  expr_ori <- expr
}

dims_use <- 1:14
resol <- 1
expr <- FindClusters(object = expr_ori, reduction.type = "pca", dims.use = dims_use, 
                     resolution = resol, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")
#PrintFindClustersParams(object = expr)
expr@meta.data$cluster <- expr@meta.data[, grep("res.", colnames(expr@meta.data), fixed = T)]
table(expr@meta.data$cluster)

# Run Non-linear dimensional reduction (tSNE)
expr <- RunTSNE(object = expr, dims.use = dims_use, nthreads = 10, do.fast = T)
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T, 
         plot.title = paste0("Dimension: 1:", max(dims_use), " Resolution: ", resol))

# check over-clustering
library("foreach")
registerDoSEQ() # avoid connection issue
expr_xxx <- ValidateClusters(expr, pc.use = dims_use, top.genes = 30, min.connectivity = 0.01, acc.cutoff = 0.85, verbose = TRUE)

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

### check batch effect
pdf(paste0(OUT, "/Seurat_batchDebug.pdf"), width = 6, height = 6, useDingbats = F)
do_batchDebug()
dev.off()
###

### check marker (expressed ratio)
marker_genes <- c("PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD14", "CD19")
do_checkMarker(x = expr_data, marker = marker_genes, cell = "")
###

### find markers for every cluster
expr.markers <- do_findMarker(expr)
table(expr.markers$cluster)
expr.markers_ftd <- expr.markers[expr.markers$power>=0.4 & expr.markers$avg_logFC>=log(2), ]
table(expr.markers_ftd$cluster)
top1 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
top5 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
###

### cell type specific marker genes
expr.markers_spec <- do_findSpecMarker(expr)
table(expr.markers_spec$cluster)
top1_spec <- expr.markers_spec %>% group_by(cluster) %>% top_n(1, avg_logFC)
top5_spec <- expr.markers_spec %>% group_by(cluster) %>% top_n(5, avg_logFC)
###

### cell-cycle gene list
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
expr <- CellCycleScoring(object = expr, s.genes = s.genes, g2m.genes = g2m.genes)
###

pdf(paste0(OUT, "/Seurat_tSNE.pdf"), width = 6, height = 6, useDingbats = F)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)
table(expr@meta.data$cluster)
library("reshape2")
cluster_distr <- acast(data = expr@meta.data[, c("bigBatch", "cluster")], formula = bigBatch ~ cluster, fun.aggregate = length, value.var = "cluster")
cluster_distr_ratio <- sweep(cluster_distr, 1, rowSums(cluster_distr), "/")
cluster_distr_ratio_melted <- melt(cluster_distr_ratio)
colnames(cluster_distr_ratio_melted) <- c("bigBatch", "cluster", "value")
ggplot(cluster_distr_ratio_melted, aes(x = factor(cluster), y = value)) + 
  geom_bar(stat = "identity", fill = scales::hue_pal()(4)[1]) + 
  facet_grid(bigBatch~.) + scale_y_continuous(expand = c(0, 0)) + xlab("Cell cluster") + ylab("Cell fraction")
if(nrow(cluster_distr) > 1) {
  pheatmap(cluster_distr, cluster_cols = F, color = colorRampPalette(c("white", "firebrick3"))(50))
}

# We also suggest exploring JoyPlot, CellPlot, and DotPlot as additional methods to view your dataset.
VlnPlot(object = expr, features.plot = top1$gene, size.title.use = 14, point.size.use = 0.1)
#VlnPlot(object = expr, features.plot = top1$gene, use.raw = TRUE, y.log = TRUE, size.title.use = 16)

for(i in top1$gene) {
  FeaturePlot(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2)
}

# DoHeatmap generates an expression heatmap for given cells and genes.
DoHeatmap(object = expr, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = expr, genes.use = c(top5$gene, "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA"), slim.col.label = TRUE, remove.key = T)
#DoHeatmap(object = expr, genes.use = c(expr.markers_ftd$gene, "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA"), slim.col.label = TRUE, remove.key = T)

top1_rmdup <- as.data.frame(top5)[! duplicated(as.data.frame(top5)$gene), ] %>% group_by(cluster) %>% top_n(1, avg_logFC)
p <- DotPlot_new(object = expr, genes.plot = top1_rmdup$gene, plot.legend = TRUE, x.lab.rot = T, rev.x = T, rev.y = T, do.plot = F, do.return = T)
p + theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))

# genes from literature
mg_literature <- read.table(paste0(OUT, "/mg_literature.txt"), header = T, sep = "\t", stringsAsFactors = F)

# cell cycle
DoHeatmap(object = expr, genes.use = c(s.genes, g2m.genes), slim.col.label = TRUE, remove.key = TRUE, cex.row = 4)
RidgePlot(object = expr, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), nCol = 2)
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = expr@meta.data$Phase)) + geom_point() + theme_bw() + 
  theme(legend.background = element_blank(), legend.box.background = element_rect(fill = alpha('white', 0.4))) + 
  theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + 
  scale_color_discrete(name = "Phase")

dev.off()

### check abnormal sub-types
pdf(file = paste0(OUT, "/Seurat_clusterDebug.pdf") , width = 5, height = 4, useDingbats = F)
do_clusterDebug()
dev.off()
###

save.image(file = paste0(OUT, "/Seurat_step2.RData"))
#load(paste0(OUT, "/Seurat_step2.RData"))

### 3. Assigning cell type identity to clusters ----
# check marker genes
data.frame(expr.markers_ftd %>% group_by(cluster) %>% top_n(2, avg_logFC), stringsAsFactors = F)
View(top10)

### GO enrichment for obtained gene markers
enriched_LS <- do_GOenrich()
topEnrich <- do.call("rbind", lapply(enriched_LS, function(x) { y <- x[[1]][1:20, c(10,1:4,9)]; y <- y[order(y$Adjusted.P.value), ] }))
View(topEnrich)
###

current.cluster.ids <- as.numeric(levels(expr@ident))
new.cluster.ids <- c("Epi-15", "Ms-E", "LUM-1", "LUM-2", "LUM-0", "Epi-13", "NC", "Ms-A", "LUM-3", "KIT", "PECAM1")
id_DF <- data.frame(old = as.character(current.cluster.ids), new = new.cluster.ids, stringsAsFactors = F)
id_DF
expr_assigned <- expr
expr_assigned@ident <- plyr::mapvalues(x = expr@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = expr_assigned, do.label = TRUE, pt.size = 1, no.legend = T)

# combine meta with identity
cellMetaDatax <- merge(cellMetaData, data.frame(ident = expr_assigned@ident, stringsAsFactors = F), by.x = 1, by.y = 0, sort = F)
cellMetaDatax$ident <- factor(cellMetaDatax$ident, levels = sort(levels(cellMetaDatax$ident)))
colnames(cellMetaDatax)[1] <- "cell"
#cellMetaDatax <- cellMetaDatax[! is.na(cellMetaDatax$ident), ]
expr.markers_ftd_labelled <- merge(x = expr.markers_ftd, y = unique(cellMetaDatax[, c("cluster", "ident")]), by = "cluster", sort = F)
expr.markers_spec_LS <- split(expr.markers_spec$gene, expr.markers_spec$cluster)
expr.markers_ftd_labelled$spec <- apply(expr.markers_ftd_labelled[, c("cluster", "gene")], 1, function(x) { 
  x[2] %in% expr.markers_spec_LS[[x[1]]]
})
#expr.markers_ftd_labelled <- expr.markers_ftd_labelled[! is.na(expr.markers_ftd_labelled$ident), ]

pdf(paste0(OUT, "/Seurat_tSNE_assigned.pdf"), width = 6, height = 6, useDingbats = F)

tmp <- expr_assigned
tmp@ident <- factor(tmp@ident, levels = sort(levels(tmp@ident)))
TSNEPlot(object = tmp, do.label = TRUE, pt.size = 1, no.legend = F)
# marker
tmp <- expr.markers_ftd_labelled %>% group_by(cluster) %>% top_n(5, avg_logFC)
tmp <- tmp[order(as.character(tmp$ident)), ]
tmp <- subset(tmp, ! is.na(ident))
DoHeatmap_new(object = expr_assigned, genes.use = tmp$gene, genes.group = tmp$ident, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, 
              group.order = sort(levels(expr_assigned@ident)), panel.spacing.y = 0, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA))) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), strip.text.y = element_blank()) + 
  theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 0, unit = "cm"))

# cell type specific marker
tmp <- subset(expr.markers_ftd_labelled, spec) %>% group_by(cluster) %>% top_n(5, avg_logFC)
tmp <- tmp[order(as.character(tmp$ident)), ]
tmp <- subset(tmp, ! is.na(ident))
DoHeatmap_new(object = expr_assigned, genes.use = tmp$gene, genes.group = tmp$ident, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, 
              group.order = sort(levels(expr_assigned@ident)), panel.spacing.y = 0, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA))) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), strip.text.y = element_blank()) + 
  theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 0, unit = "cm"))

# hierarchy marker
FeaturePlot(object = expr, features.plot = c("EPCAM", "VIM", "LMOD1", "LUM"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2)

# literature
tmp <- do.call("rbind", lapply(split(mg_literature, mg_literature$Cluster), function(x) head(x, n = 5)))
DoHeatmap_new(object = expr_assigned, genes.use = tmp$gene, genes.group = tmp$Cluster, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, 
              group.order = sort(levels(expr_assigned@ident)), panel.spacing.y = 0, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA))) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)), strip.text.y = element_blank()) + 
  theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 0, unit = "cm"))

# collagen
test_MT <- matrix(c(length(rownames(expr_data)), length(grep("^COL[0-9]", rownames(expr_data))), 
                    length(unique(expr.markers_ftd$gene)), length(unique(grep("^COL", expr.markers_ftd$gene, value = T)))), nrow = 2)
test_MT <- as.data.frame(test_MT)
test_MT <- rbind(test_MT, test_MT[1, ] - test_MT[2, ])
test_MT$non_DE <- test_MT[, 1] - test_MT[, 2]
rownames(test_MT) <- c("total", "COL", "non_COL")
colnames(test_MT) <- c("total", "DE", "non-DE")
test_MT <- test_MT[c(3,2 ), c(3, 2)]
fisher.test(test_MT[, 1:2], alternative = "greater")$p.value
sigBar <- data.frame(type=rep(c("non-DE", "DE"),each=2), ratio=c(0.011, 0.0115, 0.0115, 0.011) + 0.008)

test_DF <- as.data.frame(t(test_MT), stringsAsFactors = F)
test_DF$ratio <- test_DF[, "COL"] / rowSums(test_DF)
test_DF$type <- factor(rownames(test_DF), levels = unique(rownames(test_DF)))

ggplot(test_DF, aes(x = type, y = ratio * 100, fill = type)) + geom_bar(stat = "identity", show.legend = F) + 
  scale_y_continuous(limits = c(0, 2.2), expand = c(0, 0)) + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank()) + ylab("Ratio of Collagen genes (%)") + 
  geom_path(data = sigBar, group = 1) + geom_text(x = 1.5, y = 1.7 + 0.3, label = "***", col="red", size = 5.5) + 
  theme(plot.margin = unit(c(5.5, 60.5, 5.5, 60.5), "pt"))

tmp <- expr.markers_ftd_labelled[grep("^COL[0-9]", expr.markers_ftd_labelled$gene, value = F), ]
tmp <- tmp[order(as.character(tmp$ident)), ]
tmp <- subset(tmp, ! is.na(ident) & spec)
DoHeatmap_new(object = expr_assigned, genes.use = unique(tmp$gene), genes.group = tmp$ident[! duplicated(tmp$gene)], slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, 
              group.order = sort(levels(expr_assigned@ident)), panel.spacing.y = -0.005, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA))) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0))) + 
  theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 0, unit = "cm"))

tmp_data <- data.frame(FetchData(object = expr_assigned, vars.all = unique(tmp$gene)), check.names = F)
tmp_data$cell <- rownames(tmp_data)
tmp_data <- merge(tmp_data, cellMetaDatax[, c("cell", "ident")], by = "cell", sort = F)
tmp_data_melted <- melt(tmp_data, id.vars = c("cell", "ident"), variable.name = "gene")
tmp_data_melted$gene <- factor(tmp_data_melted$gene, levels = unique(tmp$gene))
ggplot(tmp_data_melted, aes(x = ident, y = value, fill = gene)) + geom_violin(show.legend = F) + facet_grid(gene ~ ident, scales = "free_x", switch = "y") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(strip.text.x = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.placement = "outside")

# fraction
tmp <- table(cellMetaDatax[, c("batch", "ident")], useNA = "al")
tmp_DF <- melt(tmp / rowSums(tmp))
tmp_DF <- tmp_DF[ (! is.na(tmp_DF$batch)) & (! is.na(tmp_DF$ident) ), ]
ggplot(tmp_DF, aes(x = ident, y = value * 100, fill = ident)) + 
  stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar", show.legend = F) + 
  stat_summary(fun.data=mean_cl_normal, width = 0.6, geom = "errorbar", show.legend = F) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Cell type") + ylab("Cell composition (%)") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30))

dev.off()

# check marker after assign identity
pdf(paste0(OUT, "/Seurat_tSNE_spMarker.pdf"), width = 6, height = 6, useDingbats = F)
DoHeatmap(object = expr_assigned, genes.use = c("EPCAM", "CD3D", "CD4", "CD8A", "HBA1"), 
          slim.col.label = TRUE, remove.key = TRUE, cex.row = 8) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)))

mg_DF <- data.frame(ident = expr_assigned@ident, expressed = expr@raw.data["EPCAM", names(expr@ident), drop = T] > 0, stringsAsFactors = F)
ggplot(mg_DF, aes(x = ident, fill = expressed)) + geom_bar() + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + 
  theme(legend.position = c(0.99, 0.99), legend.justification = c("right", "top")) + ylab("Cell number")

dev.off()

# 4. save and write meta table ----
save(expr, file = paste0(OUT, "/Seurat_expr.Robj"))
write.table(x = cellMetaDatax, file = paste0(OUT, "/Seurat_metaData.txt"), row.names = F, col.names = T, quote = F,sep = "\t")
write.table(x = expr.markers_ftd_labelled, file = paste0(OUT, "/Seurat_markerGenes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = expr@scale.data, file = paste0(OUT, "/UMIcount_scaled.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")

save.image(file = paste0(OUT, "/clustering.RData"))
