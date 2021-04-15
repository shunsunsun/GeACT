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
source("../../scripts/cellType_tools.r")

samplingPos <- "Epi"
OUT <- paste0("03-expression/merged/cellCluster/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/clustering.RData"))

# 1. pre-process ----
# Load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Genomes/human/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(geneID)
colnames(geneID) <- c("ensembl", "symbol")

# Load the dataset
expr_data <- read.table(file = paste0("03-expression/merged/filtering/UMIcount_filtered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F, comment.char = "")
dim(expr_data)

cell_metadata <- read.table("03-expression/merged/cellCluster/Seurat_metaData_pooled.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
dim(cell_metadata)
cell_metatable <- read.table("cell_metatable_filtered.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
dim(cell_metatable)
cell_metadata <- merge(cell_metadata, cell_metatable[, "tissue", drop = F], by = 0, sort = F)
cell_metadata$tissue <- gsub("_", " ", Hmisc::capitalize(cell_metadata$tissue))
rownames(cell_metadata) <- cell_metadata$Row.names
cell_metadata <- cell_metadata[, -1]

# add cell group
cell_metadata$clgrp <- ident2clgrp(cell_metadata$ident)
cell_metadata$ts_clgrp <- paste(cell_metadata$tissue, cell_metadata$clgrp, sep = ".")
# subset
cell_metadata_sub <- subset(cell_metadata, clgrp == "Epithelial")
expr_data_sub <- expr_data[, rownames(cell_metadata_sub)]
#

# Initialize the Seurat object with the raw (non-normalized data).
expr <- CreateSeuratObject(raw.data = expr_data_sub, min.cells = 3, min.genes = 500, project = samplingPos, names.delim = "/")
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
expr <- JackStraw(object = expr, num.pc = 50, num.replicate = 100, num.cores = 20, do.par = T)
JackStrawPlot(object = expr, PCs = 1:50)
PCElbowPlot(object = expr, num.pc = 50)

cairo_pdf(paste0(OUT, "/Seurat_dim.pdf"), width = 6, height = 10, onefile = T)
JackStrawPlot(object = expr, PCs = 1:50)
PCElbowPlot(object = expr, num.pc = 50) + theme(aspect.ratio = 1)
dev.off()

#dev.off()

save.image(file = paste0(OUT, "/Seurat_step1.RData"))
#load(paste0(OUT, "/Seurat_step1.RData"))

### 2. Cluster the cells ----
if(! exists("expr_ori")) {
  print("Create copy for original expr")
  expr_ori <- expr
}

dims_use <- 1:40  # 1:25
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
all(rownames(expr@meta.data) == rownames(cell_metadata_sub))
expr@meta.data$cluster <- cell_metadata_sub$tissue
all(names(expr@ident) == rownames(cell_metadata_sub))
expr@ident <- factor(cell_metadata_sub$tissue)
names(expr@ident) <- rownames(cell_metadata_sub)
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
cellMetaData$cluster <- factor(cellMetaData$cluster, levels = sort(unique(cellMetaData$cluster)))
cellStat <- read.table(paste0("03-expression/merged/filtering/", ".", "/filtering_cells.txt"), header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
cellMetaData <- merge(cellMetaData, cellStat[, -c(4, 5)], by.x = 1, by.y = 0, sort = F)
colnames(cellMetaData)[1] <- "cell"
# add tissue/samplingPos/ident
cell_metatable <- read.table("cell_metatable.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cell_metatable)
cellMetaData <- merge(cellMetaData, cell_metatable[, c("cell", "tissue", "samplingPos", "ident")], by = "cell", sort = F)
cellMetaData$tissue <- gsub("_", " ", Hmisc::capitalize(cellMetaData$tissue))
cellMetaData$tissue <- factor(cellMetaData$tissue, levels = sort(unique(cellMetaData$tissue)))
cellMetaData$ident_old <- cellMetaData$ident
cellMetaData$ident <- cellMetaData$tissue

### check batch effect
# cairo_pdf(paste0(OUT, "/Seurat_batchDebug.pdf"), width = 6, height = 6, onefile = T)
# do_batchDebug()
# dev.off()
###

### check marker (expressed ratio)
# marker_genes <- c("PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD14", "CD19")
# do_checkMarker(x = expr_data, marker = marker_genes, cell = "")
###

## find markers for every cluster
expr.markers_pn <- do_findMarker(expr, only.pos = F, ncpu = 12)
expr.markers <- subset(expr.markers_pn, avg_logFC > 0)
table(expr.markers$cluster)
expr.markers_ftd <- expr.markers[expr.markers$power>=0.4 & expr.markers$avg_logFC>=log(2), ]
table(expr.markers_ftd$cluster)
top1 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
top3 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(3, avg_logFC)
top5 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
##

## cell type specific marker genes
expr.markers_spec <- do_findSpecMarker(expr, ncpu = 12)
table(expr.markers_spec$cluster)
top1_spec <- expr.markers_spec %>% group_by(cluster) %>% top_n(1, avg_logFC)
top5_spec <- expr.markers_spec %>% group_by(cluster) %>% top_n(5, avg_logFC)
##

### cell-cycle gene list
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
expr <- CellCycleScoring(object = expr, s.genes = s.genes, g2m.genes = g2m.genes)
###

# the ident of each cluster
cluster_ident_fst5 <- t(sapply(stringr::str_sort(unique(expr@meta.data$cluster), numeric = T), function(x) {
  y0 <- sort(table(subset(cell_metatable, cell %in% rownames(expr@meta.data)[expr@meta.data$cluster == x], "ident_old", drop = T)), decreasing = T)
  y <- paste0(names(y0), "(", y0, ")")[1:5]
  return(y)
}))
cluster_ident_fst5 <- as.data.frame(cluster_ident_fst5)
View(cluster_ident_fst5)

pdf(paste0(OUT, "/Seurat_tSNE.pdf"), width = 6, height = 6, useDingbats = F)

TSNEPlot(object = expr, pt.size = 1, do.label = T, no.legend = T)

# color by tissue
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = tissue)) + geom_point(size = 1) +
  #theme(legend.position = "bottom", legend.justification = c("center")) + 
  theme(aspect.ratio = 1) + 
  #theme(legend.box.margin = margin(t = -15), plot.margin = margin(t = 10, l = -20, r = 20)) + 
  #theme(legend.key.width = unit(0.2, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") #+ 
  #guides(color = guide_legend(ncol = 4, override.aes = list(size = 3.6))) + 
  #scale_color_manual(values = c(brewer.pal(name = "Dark2", n = 8), brewer.pal(name = "Paired", n = 9)))

dev.off()

### check abnormal sub-types
# pdf(file = paste0(OUT, "/Seurat_clusterDebug.pdf") , width = 5, height = 4, useDingbats = F)
# do_clusterDebug()
# dev.off()
###

save.image(file = paste0(OUT, "/Seurat_step2.RData"))
#load(paste0(OUT, "/Seurat_step2.RData"))

### 3. Assigning cell type identity to clusters ----
# check marker genes
data.frame(expr.markers_ftd %>% group_by(cluster) %>% top_n(2, avg_logFC), stringsAsFactors = F)
View(top10)

### GO enrichment for obtained gene markers
enriched_LS <- do_GOenrich(expr.markers_ftd, ncpu = 12)
topEnrich <- do.call("rbind", lapply(enriched_LS, function(x) { y <- x[1:20, c(9,1:8)] }))
View(topEnrich)
###

current.cluster.ids <- levels(expr@ident)
new.cluster.ids <- levels(expr@ident)

id_DF <- data.frame(old = as.character(current.cluster.ids), new = new.cluster.ids, stringsAsFactors = F)
id_DF
expr_assigned <- expr
# expr_assigned@ident <- plyr::mapvalues(x = expr@ident, from = current.cluster.ids, to = new.cluster.ids)
# # change sort
cell_type_sorted <- new.cluster.ids[c(2,10,9,4,5,8,1,6,3,7)]
length(cell_type_sorted)
expr_assigned@ident <- factor(expr_assigned@ident, levels = cell_type_sorted)

TSNEPlot(object = expr_assigned, do.label = TRUE, pt.size = 1, no.legend = T)

# combine meta with identity
cellMetaDatax <- cellMetaData
cellMetaDatax$ident <- factor(cellMetaDatax$ident, levels = cell_type_sorted)
#colnames(cellMetaDatax)[1] <- "cell"
#cellMetaDatax <- cellMetaDatax[! is.na(cellMetaDatax$ident), ]
expr.markers_ftd_labelled <- merge(x = expr.markers_ftd, y = unique(cellMetaDatax[, c("cluster", "ident")]), by = "cluster", sort = F)
expr.markers_spec_LS <- split(expr.markers_spec$gene, expr.markers_spec$cluster)
expr.markers_ftd_labelled$spec <- apply(expr.markers_ftd_labelled[, c("cluster", "gene")], 1, function(x) { 
  x[2] %in% expr.markers_spec_LS[[x[1]]]
})
#expr.markers_ftd_labelled <- expr.markers_ftd_labelled[! is.na(expr.markers_ftd_labelled$ident), ]

# color set
load("03-expression/merged/cellCluster/ct_color.RData")
ct_color
ct_color <- ct_color[levels(expr_assigned@ident)]

# # stat by ident
# # general correlation
# all(colnames(expr@raw.data) == cellMetaDatax$cell)
# tmp_LS <- split(as.data.frame(t(expr@raw.data)), cellMetaDatax$ident)
# 
# # average cor
# cl <- makeCluster(min(length(levels(expr@ident)), 12))
# clusterExport(cl = cl, varlist = "tmp_LS")
# cor_avg <- parSapply(cl, tmp_LS, function(x) {
#   x1 <- t(x)
#   y <- sapply(tmp_LS, function(k) {
#     x2 <- t(k)
#     ky <- mean(cor(x1, x2))
#   })
# })
# stopCluster(cl); rm(cl)
# 
# # average profile
# tmp <- sapply(tmp_LS, function(x) { y <- colSums(x); y <- y / sum(y); return(y) })
# 
# tmp_cor <- cor_avg
# tmp_cor[lower.tri(tmp_cor)] <- NA
# #corrplot(tmp_cor, type="upper", method = "color", col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
# tmp_cor_melted <- na.omit(melt(tmp_cor))
# tmp_cor_melted$Var2 <- factor(tmp_cor_melted$Var2, levels = rev(levels(tmp_cor_melted$Var2)))
# p1 <- ggplot(tmp_cor_melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
#   scale_fill_gradientn(name = "Correlation", colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
#   #geom_text(data = subset(tmp_cor_melted, Var1 == Var2), aes(label = Var1), angle = 45, hjust = 0, nudge_x = -1.5, nudge_y = -1.5) + 
#   theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
#   theme(aspect.ratio = 1, legend.position = c(0.8, 0.88))
# 
# tmp_cor <- cor(tmp)
# tmp_cor[lower.tri(tmp_cor)] <- NA
# #corrplot(tmp_cor, type="upper", method = "color", col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
# tmp_cor_melted <- na.omit(melt(tmp_cor))
# tmp_cor_melted$Var2 <- factor(tmp_cor_melted$Var2, levels = rev(levels(tmp_cor_melted$Var2)))
# 
# p2 <- ggplot(tmp_cor_melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
#   scale_fill_gradientn(name = "Correlation", colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
#   #geom_text(data = subset(tmp_cor_melted, Var1 == Var2), aes(label = Var1), angle = 45, hjust = 0, nudge_x = -1.5, nudge_y = -1.5) + 
#   theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
#   theme(aspect.ratio = 1, legend.position = c(0.8, 0.88))
# 
# pdf(file = paste0(OUT, "/", samplingPos, "/Seurat_cellType_heatmap.pdf"), width = 8, height = 8)
# print(p1)
# print(p2)
# dev.off()

# plot GO enrich
enrich_res <- do.call("rbind", enriched_LS)
enrich_res$Term <- Hmisc::capitalize(gsub(" \\(.*\\)", "", enrich_res$Term))
enrich_res$Adjusted.P.value <- enrich_res$q_value
#enrich_res <- enrich_res[order(enrich_res$Adjusted.P.value), ]
enrich_res <- subset(enrich_res, Adjusted.P.value < 0.05)
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res, enrich_res$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, id_DF, by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# add dummy
nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new), 
                                                 Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
enrich_res_DF$new <- factor(enrich_res_DF$new, levels = cell_type_sorted)
enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", "/enrichment_all.pdf"), width = 8, height = 8)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) + 
  geom_bar(stat = "identity", width = 0.9, show.legend = F) + 
  facet_grid(new ~ ., scales ="free", drop = F) + 
  coord_flip() + 
  ylab(expression(paste(-log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.01)) + 
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") + 
  theme(aspect.ratio = 0.4) +
  scale_fill_manual(values = ct_color, drop = F)
gp
dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_assigned.pdf"), width = 6, height = 6, useDingbats = F) # 6*6

tmp <- expr_assigned
#tmp@ident <- factor(tmp@ident, levels = sort(levels(tmp@ident)))
ident_labels <- paste(1:length(levels(tmp@ident)), levels(tmp@ident), sep = ": ")
p <- TSNEPlot(object = tmp, do.label = TRUE, pt.size = 1, label.size = 3, no.legend = F, do.return = T) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  scale_color_manual(values = ct_color) + 
  theme(aspect.ratio = 1, panel.border = element_blank()) + 
  #guides(color = guide_legend(ncol = 1)) + 
  coord_cartesian(clip = "off") + 
  xlab("tSNE-1") + ylab("tSNE-2")
###
p$layers[[3]]$data[7, 2:3] <- list(3, -33)
p$layers[[3]]$data[9, 2] <- -30
p$layers[[3]]$data[4, 3] <- -1
p$layers[[3]]$data[10, 3] <- -5
p$layers[[3]]$data[2, 2:3] <- list(12, -7)
###
p

# # color by seq ID
# cellSeqID <- read.table("01-cleandata/merged/cleanFqStat.txt", header = F, sep = "\t", stringsAsFactors = F)[, c(5, 13)]
# dim(cellSeqID)
# colnames(cellSeqID) <- c("cell", "seqid")
# tmp@meta.data$seqid <- cellSeqID[match(rownames(tmp@meta.data), cellSeqID$cell), "seqid"]
# 
# DimPlot(object = tmp, reduction.use = "tsne", group.by = "seqid", do.label = F, pt.size = 1, no.legend = F, do.return = T) + 
#   theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank(), legend.background = element_blank()) + 
#   theme(legend.position = c(0.135, 0.10)) + 
#   xlab("tSNE-1") + ylab("tSNE-2")

# color by bigBatch
# if(length(unique(expr@meta.data$bigBatch)) > 1) {
#   ### change levels
#   tmp@meta.data$bigBatch <- factor(tmp@meta.data$bigBatch, levels = c("X"))
#   levels(tmp@meta.data$bigBatch)
#   levels(tmp@meta.data$bigBatch) <- c("X")
#   levels(tmp@meta.data$bigBatch)
#   ###
#   DimPlot(object = tmp, reduction.use = "tsne", group.by = "bigBatch", do.label = F, pt.size = 1, no.legend = F, do.return = T) + 
#     theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank(), legend.background = element_blank()) + 
#     theme(legend.position = c(0.125, 0.10)) + 
#     xlab("tSNE-1") + ylab("tSNE-2")
# }

# hierarchy marker
FeaturePlot_new(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "HBG1"), new.title = c("EPCAM (CD326)", "VIM", "PTPRC (CD45)", "HBG1"), 
                cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 1, no.axes = T)

# dot plot
genes_list <- c("KRT15", "CLDN18", "MUC13", "LGALS4", "MYB", "CLPS", "ANXA1", "CPM", "PODXL", "CCDC170")
p <- DotPlot_new(object = expr_assigned, genes.plot = genes_list, 
                 plot.legend = TRUE, x.lab.rot = T, rev.x = T, rev.y = T, do.plot = F, do.return = T, breaks.use = c(-1, 0, 1, 2), limits.use = c(-1.5, 2.5)) + 
  theme(aspect.ratio = 0.525) + 
  theme(legend.position = "bottom", legend.justification = "center") + 
  theme(legend.box.margin = margin(t = -15, r = 20)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5, linetype = 1))
# p$layers <- c(geom_rect(xmin = 2-0.5, xmax = 2+0.5, ymin = 1-0.5, ymax = 4+0.5, fill = "skyblue"), p$layers)
# p$layers <- c(geom_rect(xmin = 3-0.5, xmax = 3+0.5, ymin = 1-0.5, ymax = 4+0.5, fill = "skyblue"), p$layers)
p + guides(fill = guide_colorbar(title = "Z-score", order = 1), size = guide_legend(title = "Percentage", label.vjust = 1.95, label.position = "bottom", override.aes = list(colour = "black")))

# proliferative score (cell cycle)
tmp <- merge(expr_assigned@meta.data, id_DF, by.x = "cluster", by.y = "old", sort = F)
tmp$new <- factor(tmp$new, levels = cell_type_sorted)
ggplot(tmp, aes(x = new, y = S.Score + G2M.Score, fill = new)) + geom_boxplot(show.legend = F) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Proliferative score") + 
  scale_fill_manual(values = ct_color, drop = F)

# TF
tf_list <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/AnimalTFDB/human_TF_list.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]
tmp <- subset(expr.markers_ftd_labelled, gene %in% tf_list, drop = F)
tmp <- tmp[order(as.character(tmp$ident)), ]
tmp <- tmp[! duplicated(tmp$gene), ]
tmp <- do.call("rbind", split(tmp, tmp$ident))
tmp <- do.call("rbind", lapply(split(tmp, tmp$ident), function(x) { head(x, 1) }))

# violin plot
tmp_data <- data.frame(FetchData(object = expr_assigned, vars.all = unique(tmp$gene)), check.names = F)
tmp_data$cell <- rownames(tmp_data)
tmp_data <- merge(tmp_data, cellMetaDatax[, c("cell", "ident")], by = "cell", sort = F)
tmp_data_melted <- melt(tmp_data, id.vars = c("cell", "ident"), variable.name = "gene")
tmp_data_melted$ident <- factor(tmp_data_melted$ident, levels = rev(levels(tmp_data_melted$ident)))
tmp_data_melted$gene <- factor(tmp_data_melted$gene, levels = unique(tmp$gene))
p <- ggplot(tmp_data_melted, aes(x = ident, y = value, fill = ident, color = ident)) + geom_violin(show.legend = F, scale = "width") + 
  theme(aspect.ratio = 6.25) + 
  facet_grid(. ~ gene, scales = "free_x", switch = "x") + coord_flip() + 
  scale_fill_manual(values = ct_color, drop = F) + 
  scale_color_manual(values = ct_color, drop = F) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  theme(panel.grid = element_blank()) + 
  theme(strip.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = -10, r = 10, b = 15)), strip.background = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linetype = 1, size = 0.75))
pg <- ggplotGrob(p)
for(i in which(grepl("strip-b", pg$layout$name))){ pg$grobs[[i]]$layout$clip <- "off" }
grid.newpage(); grid::grid.draw(pg)

# collagen
test_MT <- matrix(c(length(rownames(expr_data)), length(grep("^COL[0-9]", rownames(expr_data))), 
                    length(unique(expr.markers_ftd$gene)), length(unique(grep("^COL[0-9]", expr.markers_ftd$gene, value = T)))), nrow = 2)
test_MT <- as.data.frame(test_MT)
test_MT <- rbind(test_MT, test_MT[1, ] - test_MT[2, ])
test_MT$non_DE <- test_MT[, 1] - test_MT[, 2]
rownames(test_MT) <- c("total", "COL", "non_COL")
colnames(test_MT) <- c("total", "DE", "non-DE")
test_MT <- test_MT[c(3,2 ), c(3, 2)]
pv <- fisher.test(test_MT[, 1:2], alternative = "greater")$p.value
print(pv)
if(pv < 0.001) {
  sigSym <- "***"
} else if(pv < 0.01) {
  sigSym <- "**"
} else if(pv < 0.05) {
  sigSym <- "*"
} else {
  sigSym <- "N.S."
}
sigBar <- data.frame(type=rep(c("non-DE", "DE"),each=2), ratio=c(0.011, 0.0115, 0.0115, 0.011) + 0.008)

test_DF <- as.data.frame(t(test_MT), stringsAsFactors = F)
test_DF$ratio <- test_DF[, "COL"] / rowSums(test_DF)
test_DF$type <- factor(rownames(test_DF), levels = unique(rownames(test_DF)))

ggplot(test_DF, aes(x = type, y = ratio * 100, fill = type)) + geom_bar(stat = "identity", show.legend = F) + 
  scale_y_continuous(limits = c(0, 2.2), expand = c(0, 0)) + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank()) + ylab("Ratio of Collagen genes (%)") + 
  geom_path(data = sigBar, group = 1) + geom_text(x = 1.5, y = 1.7 + 0.3, label = sigSym, col="red", size = 5.5) + 
  theme(plot.margin = unit(c(5.5, 60.5, 5.5, 60.5), "pt"))

tmp <- expr.markers_ftd_labelled[grep("^COL[0-9]", expr.markers_ftd_labelled$gene, value = F), ]
tmp <- tmp[order(as.character(tmp$ident)), ]
#tmp <- subset(tmp, ! is.na(ident) & spec)
tmp <- tmp[! duplicated(tmp$gene), ]
tmp <- do.call("rbind", split(tmp, tmp$ident))

# violin plot
if(nrow(tmp) > 0) {
  tmp_data <- data.frame(FetchData(object = expr_assigned, vars.all = unique(tmp$gene)), check.names = F)
  tmp_data$cell <- rownames(tmp_data)
  tmp_data <- merge(tmp_data, cellMetaDatax[, c("cell", "ident")], by = "cell", sort = F)
  tmp_data_melted <- melt(tmp_data, id.vars = c("cell", "ident"), variable.name = "gene")
  tmp_data_melted$ident <- factor(tmp_data_melted$ident, levels = rev(levels(tmp_data_melted$ident)))
  tmp_data_melted$gene <- factor(tmp_data_melted$gene, levels = unique(tmp$gene))
  ggplot(tmp_data_melted, aes(x = ident, y = value, fill = ident, color = ident)) + geom_violin(show.legend = F, scale = "width") + 
    theme(aspect.ratio = 6.25) + 
    facet_grid(. ~ gene, scales = "free_x", switch = "x") + coord_flip() + 
    scale_fill_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident)) + 1)[-1]), drop = F) + 
    scale_color_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident)) + 1)[-1]), drop = F) + 
    theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
    theme(panel.grid = element_blank()) + 
    theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5)), strip.background = element_blank()) + 
    theme(panel.border = element_rect(colour = "black", fill = NA, linetype = 1, size = 0.75))
}

# fraction
# tmp <- table(cellMetaDatax[, c("batch", "ident")], useNA = "al")
# tmp_DF <- melt(tmp / rowSums(tmp))
# tmp_DF <- tmp_DF[ (! is.na(tmp_DF$batch)) & (! is.na(tmp_DF$ident) ), ]
# ggplot(tmp_DF, aes(x = ident, y = value * 100, fill = ident)) + 
#   stat_summary(fun.y=mean,position=position_dodge(width=0.95), geom="bar", show.legend = F) + 
#   stat_summary(fun.data=mean_cl_normal, width = 0.6, geom = "errorbar", show.legend = F) + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Cell type") + ylab("Cell type composition (%)") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp_DF$value) * 100 * 1.05))

# # fraction (pooled)
# tmp <- table(factor(subset(cellMetaDatax, ! ident %in% c("FGC", "FGC/Granulosa"), "ident", drop = T)), useNA = "al")
# tmp_DF <- melt(tmp / sum(tmp))
# colnames(tmp_DF)[1] <- "ident"
# tmp_DF <- tmp_DF[ ! is.na(tmp_DF$ident), ]
# ggplot(tmp_DF, aes(x = ident, y = value, fill = ident)) + 
#   geom_bar(stat = "identity", show.legend = F) + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Cell type") + ylab("Cell type composition") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp_DF$value) * 1.05)) + 
#   scale_fill_manual(values = scales::hue_pal()(length(levels(tmp_DF$ident)) + 1)[-1], drop = F)
# 
# # fraction (pooled - pie)
# tmp_x <- rep(1.6, nrow(tmp_DF))
# tmp_x[8:11] <- seq(1.7, 2, length.out = 4)
# ggplot(tmp_DF, aes(x = factor(1), y = value, fill = ident)) + 
#   geom_bar(stat = "identity", width = 1, show.legend = F, color = "white") + 
#   # geom_path(data = data.frame(x = c(1.55, 1.55), y = 1 - c(cumsum(tmp_DF$value)[grep("^Fibro-", tmp_DF$ident)[1] - 1], cumsum(tmp_DF$value)[rev(grep("^Fibro-", tmp_DF$ident))[1]])),
#   #           aes(x = x, y = y), color = "#00BF74", size = 1.5, inherit.aes = F) +
#   geom_path(data = data.frame(x = c(1.55, 1.55), y = 1 - c(cumsum(tmp_DF$value)[1], cumsum(tmp_DF$value)[5])),
#             aes(x = x, y = y), color = "mediumpurple1", size = 1.5, inherit.aes = F) +
#   coord_polar(theta = "y", direction = -1) + 
#   #geom_text(aes(label = ident, color = ident, x = tmp_x, y = 1 - (cumsum(tmp_DF$value) - tmp_DF$value/2), angle = -(cumsum(tmp_DF$value) - tmp_DF$value/2) * 360), show.legend = F) + 
#   theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
#   scale_fill_manual(values = scales::hue_pal()(length(levels(tmp_DF$ident)) + 1)[-1], drop = F)
# rm(tmp_x)

# fraction in different pos
# tmp <- table(cellMetaDatax[, c("batch", "ident", "bigBatch")], useNA = "al")
# tmp_DF <- melt(tmp / rowSums(tmp))
# tmp_DF <- tmp_DF[tmp_DF$bigBatch == gsub("-.*", "", tmp_DF$batch), ]
# tmp_DF <- tmp_DF[ (! is.na(tmp_DF$batch)) & (! is.na(tmp_DF$ident) ), ]
# # sort by real position
# tmp_DF$bigBatch <- factor(tmp_DF$bigBatch, levels = c("FZY_C", "FZW_C"))
# ### change levels
# levels(tmp_DF$bigBatch)
# levels(tmp_DF$bigBatch) <- c("Central", "Surrounding")
# levels(tmp_DF$bigBatch)
# ###
# ggplot(tmp_DF, aes(x = ident, y = value * 100, fill = bigBatch)) + 
#   stat_summary(fun.y=mean, position=position_dodge(width=0.8), width = 0.75, geom="bar", show.legend = T) + 
#   stat_summary(fun.data=mean_se, position=position_dodge(width=0.8), width = 0.6, geom = "errorbar", show.legend = F) + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Cell type") + ylab("Cell type composition (%)") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + 
#   theme(legend.title = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top"))

# # fraction in different pos (bar with connected lines)
# cellMetaDatas <- cellMetaDatax
# cellMetaDatas$bigBatch <- gsub("_[A-Z]$", "", as.character(cellMetaDatas$bigBatch))
# tmp <- table(subset(cellMetaDatas, T, c("bigBatch", "ident")))
# tmp_DF <- melt(tmp / rowSums(tmp))
# # sort by real position
# tmp_DF$bigBatch <- factor(tmp_DF$bigBatch, levels = c("XX"))
# 
# ### change levels
# levels(tmp_DF$bigBatch)
# levels(tmp_DF$bigBatch) <- c("mixture")
# levels(tmp_DF$bigBatch)
# 
# tmp_DF_sub <- subset(tmp_DF, bigBatch == rev(levels(tmp_DF$bigBatch))[1])
# 
# ggplot(tmp_DF, aes(x = bigBatch, y = value, fill = ident)) + geom_bar(stat = "identity", width = 0.8, color = "white") + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank()) + 
#   scale_y_continuous(expand = c(0, 0)) + 
#   # geom_segment(mapping = aes(x = 1+0.8/2, xend = 2-0.8/2, y = 1, yend = 1), size = 0.3) + 
#   # geom_segment(data = dcast(tmp_DF, ident ~ bigBatch), mapping = aes(x = 1+0.8/2, xend = 2-0.8/2, y = (1-cumsum(IV)), yend = (1-cumsum(V))), size = 0.3) + 
#   # geom_segment(mapping = aes(x = 2+0.8/2, xend = 3-0.8/2, y = 1, yend = 1), size = 0.3) + 
#   # geom_segment(data = dcast(tmp_DF, ident ~ bigBatch), mapping = aes(x = 2+0.8/2, xend = 3-0.8/2, y = (1-cumsum(V)), yend = (1-cumsum(VIII))), size = 0.3) + 
#   # geom_path(data = data.frame(x = c(4.5, 4.5), y = 1 - c(cumsum(tmp_DF_sub$value)[grep("^Fibro-", tmp_DF_sub$ident)[1] - 1], cumsum(tmp_DF_sub$value)[rev(grep("^Fibro-", tmp_DF_sub$ident))[1]])),
#   #           aes(x = x, y = y), color = "#00BF74", size = 1.5, inherit.aes = F) +
#   geom_path(data = data.frame(x = c(1.5, 1.5), y = 1 - c(cumsum(tmp_DF_sub$value)[1], cumsum(tmp_DF_sub$value)[5])),
#             aes(x = x, y = y), color = "mediumpurple1", size = 1.5, inherit.aes = F) + 
#   # guides(fill = guide_legend(ncol = 1)) + 
#   ylab("Cell type composition") + 
#   scale_fill_manual(values = scales::hue_pal()(length(levels(tmp_DF$ident)) + 1)[-1], drop = F)

dev.off()

# 4. save and write meta table ----
#save(expr, file = paste0(OUT, "/Seurat_expr.Robj"))
write.table(x = cellMetaDatax, file = paste0(OUT, "/Seurat_metaData.txt"), row.names = F, col.names = T, quote = F,sep = "\t")
write.table(x = expr.markers_pn, file = paste0(OUT, "/Seurat_markerGenes_pn.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = expr.markers_ftd_labelled, file = paste0(OUT, "/Seurat_markerGenes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
data.table::fwrite(x = as.data.frame(expr@scale.data), file = paste0(OUT, "/UMIcount_scaled.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)

save.image(file = paste0(OUT, "/clustering.RData"))
