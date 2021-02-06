# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

suppressMessages(library("arrow"))
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
library("Matrix")
library("pheatmap")
library("reshape2")
library("grid")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
suppressMessages(library("topGO"))
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
expr_data <- read_feather(file = file.path("03-expression/merged/filtering", samplingPos, "UMIcount_filtered.feather"))
expr_data <- as.data.frame(expr_data)
expr_data_gene <- read.table(file = file.path("03-expression/merged/filtering", samplingPos, "UMIcount_filtered.gene"), header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_data) <- expr_data_gene$V1

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
expr <- JackStrawPlot(object = expr, PCs = 1:100)
PCElbowPlot(object = expr, num.pc = 100)

dim_limit <- min(which(expr@dr$pca@jackstraw@overall.p.values[, 2] > 0.001)) - 1
cat("Suggested max. dimension:", dim_limit, "\n")

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

dims_use <- 1:dim_limit  # 1:25
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
expr <- RunUMAP(object = expr, reduction.use = "pca", dims.use = dims_use, min_dist = 0.75)

DimPlot(object = expr, reduction.use = "umap", no.legend = T, do.return = TRUE, vector.friendly = TRUE, pt.size = 3) + 
  ggtitle(paste0("Dimension: 1:", max(dims_use), " Resolution: ", resol)) + theme(plot.title = element_text(hjust = 0.5))

# add info to meta
cellMetaData <- merge(expr@meta.data, cbind(expr@dr$tsne@cell.embeddings, expr@dr$umap@cell.embeddings), by = 0, sort = F)
cellMetaData$batch <- factor(cellMetaData$batch, levels = unique(cellMetaData$batch))
cellMetaData$cluster <- factor(cellMetaData$cluster, levels = sort(as.numeric(unique(cellMetaData$cluster))))
cellStat <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/filtering_cells.txt"), header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
cellMetaData <- merge(cellMetaData, cellStat[, -c(4, 5)], by.x = 1, by.y = 0, sort = F)
colnames(cellMetaData)[1] <- "cell"
# add tissue/samplingPos/ident/group
cell_metatable <- read.table("cell_metatable_filtered_plus.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cell_metatable)
cellMetaData <- merge(cellMetaData, cell_metatable[, c("cell", "tissue", "samplingPos", "ident", "group")], by = "cell", sort = F)
cellMetaData$tissue <- gsub("_", " ", Hmisc::capitalize(cellMetaData$tissue))
cellMetaData$tissue <- factor(cellMetaData$tissue, levels = sort(unique(cellMetaData$tissue)))
cellMetaData$group_new <- cellMetaData$group
group_shown <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "B", 
                 "DC/Macrophage", "NKT", "T", "Glial", "FGC", 
                 "Granulosa", "Sertoli", "Erythrocyte", "Other")
cellMetaData$group_new[! cellMetaData$group_new %in% group_shown] <- "Other"
cellMetaData$group_new <- factor(cellMetaData$group_new, levels = group_shown)

write.table(x = cellMetaData[, colnames(cellMetaData) != "orig.ident"], file = "cell_metatable_merged.txt", row.names = F, col.names = T, quote = F, sep = "\t")

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

# set color by tissue
ct_color <- c(brewer.pal(name = "Dark2", n = 8), brewer.pal(name = "Paired", n = 9))
names(ct_color) <- levels(cellMetaData$tissue)
ct_color["Kidney"] <- "#F2E100"
ct_color["Large intestine"] <- "tomato"
ct_color["Lung"] <- "#419FDE"
save(ct_color, file = paste0(OUT, "/ct_color.RData"))

# set color by cell group
cg_color <- c("tomato", "orange", "#FB9A99", "#00BF74", "#8258FA", 
              "#D0A9F5", "purple2", "#D358F7", "#419FDE", "#F2E100", 
              ct_color["Ovary"], ct_color["Testis"], "firebrick3", "grey80")
names(cg_color) <- group_shown
save(cg_color, file = paste0(OUT, "/cg_color.RData"))

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

pdf(paste0(OUT, "/Seurat_tSNE_legend_1col.pdf"), width = 9, height = 6, useDingbats = F)

# color by tissue
ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = tissue)) + geom_point(size = 0.2) +
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.8, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color)

dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_and_legend.pdf"), width = 6, height = 5, useDingbats = F)

# color by tissue
gp <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = tissue)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.725, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color) + 
  ggtitle("Transcriptome") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

# color by cell group
gp <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = group_new)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.8, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = cg_color) + 
  ggtitle("Transcriptome") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

dev.off()

pdf(paste0(OUT, "/Seurat_UMAP_and_legend.pdf"), width = 6, height = 5, useDingbats = F)

# color by tissue
gp <- ggplot(cellMetaData, aes(x = UMAP1, y = UMAP2, color = tissue)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.725, "cm")) + 
  labs(color = NULL) + xlab("UMAP1") + ylab("UMAP2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color) + 
  ggtitle("Transcriptome") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

# color by cell group
gp <- ggplot(cellMetaData, aes(x = UMAP1, y = UMAP2, color = group_new)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.8, "cm")) + 
  labs(color = NULL) + xlab("UMAP1") + ylab("UMAP2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = cg_color) + 
  ggtitle("Transcriptome") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

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

png(paste0(OUT, "/Seurat_tSNE_simple.png"), bg = "transparent", width = 3000, height = 3000, res = 600)

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

pdf(paste0(OUT, "/Seurat_tSNE_spMarker.pdf"), width = 4, height = 5, useDingbats = F)

# FeaturePlot
mk_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1", "HBG1")
mk_labels <- c("EPCAM (CD326)", "PECAM1 (CD31)", "PTPRC (CD45)", "COL1A1", "HBG1")
for(i in seq_along(mk_genes)) {
  p <- FeaturePlot_new(object = expr, features.plot = mk_genes[i], new.title = mk_labels[i], cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 1, no.axes = T, do.plot = F, do.return = T)[[1]]
  p <- p + theme(aspect.ratio = 1, plot.title = element_text(size = 20), plot.margin = margin(5,2,-10,2))
  print(p)
}

dev.off()

### signature genes
ident.ori <- expr@ident
expr_assigned <- expr

cell_metatable <- cellMetaData
rownames(cell_metatable) <- cell_metatable$cell
cell_metatable$ts_ident <- Hmisc::capitalize(paste(cell_metatable$tissue, cell_metatable$ident, sep = "."))
# sort
cell_metatable <- cell_metatable[order(cell_metatable$tissue), ]
cell_metatable_LS <- split(cell_metatable, cell_metatable$group_new)
names(cell_metatable_LS) <- NULL
cell_metatable_sorted <- do.call("rbind", cell_metatable_LS)
#
cell_type_ti <- cell_metatable[names(ident.ori), "ts_ident"]
names(cell_type_ti) <- names(ident.ori)
cell_type_ti <- factor(cell_type_ti, levels = unique(cell_metatable_sorted$ts_ident))
levels(cell_type_ti) <- 1:length(levels(cell_type_ti))
expr_assigned@ident <- cell_type_ti

signatureGenes <- read.table("../../pooled_data/All/cellType_signatureGenes.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(signatureGenes)
signatureGenes$ts_ident <- Hmisc::capitalize(paste(signatureGenes$tissue, signatureGenes$ident, sep = "."))
signatureGenes <- merge(x = signatureGenes, y = unique(cell_metatable[, c("ts_ident", "group_new")]), by = "ts_ident", sort = F)
signatureGenes <- signatureGenes[, c(2:13,1,14)]
###
signatureGenes <- subset(signatureGenes, ! group_new %in% "Other")
###
signatureGenes_LS <- split(signatureGenes, signatureGenes$ts_ident)
signatureGenes_LS <- signatureGenes_LS[unique(cell_metatable_sorted$ts_ident)]
signatureGenes_LS <- lapply(signatureGenes_LS, function(x) { head(x, 2) })
names(signatureGenes_LS) <- NULL
signatureGenes_sorted <- do.call("rbind", signatureGenes_LS)
signatureGenes_unique <- signatureGenes_sorted[! duplicated(signatureGenes_sorted$gene), ]
gene_hlt <- unlist(lapply(split(signatureGenes_unique$gene, signatureGenes_unique$group_new), function(x) { head(x, 1) }))

data_use <- DoHeatmap_object(object = expr_assigned, genes.use = signatureGenes_unique$gene, slim.col.label = TRUE, remove.key = F, rotate.key = T, do.plot = F,
                             group.spacing = 0,
                             cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA), max.cells.per.ident = 3, random.seed = 1))
data_use_MT <- acast(data_use[, 1:3], formula = gene ~ cell, value.var = "expression")
data_use_MT <- data_use_MT[signatureGenes_unique$gene, ]
cell_metatable_sel <- merge(x = data.frame(cell = colnames(data_use_MT), stringsAsFactors = F), y = cellMetaData[, c("cell", "tissue", "group_new")], by = "cell", sort = F)
###
cell_metatable_sel <- subset(cell_metatable_sel, group_new != "Other")
data_use_MT <- data_use_MT[, colnames(data_use_MT) %in% cell_metatable_sel$cell]
###

pdf(paste0(OUT, "/Seurat_tSNE_assigned_heatmap.pdf"), width = 10, height = 6, useDingbats = F)

# DoHeatmap(object = expr_assigned, genes.use = signatureGenes_unique$gene, slim.col.label = TRUE, remove.key = F, rotate.key = T, do.plot = F, 
#           group.spacing = 0, 
#           cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA), max.cells.per.ident = 3, random.seed = 1)) + 
#   theme(axis.text.y = element_blank(), strip.text.x = element_blank()) + 
#   theme(legend.position = "bottom", legend.justification = "center", legend.box.spacing = unit(0.25, units = "cm")) + 
#   guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(4, units = "cm")))


library("ComplexHeatmap")
col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#4575B4", "white", "#D73027"))

ht_opt(legend_title_gp = gpar(fontsize = 12), legend_labels_gp = gpar(fontsize = 12))
ht1 <- Heatmap(matrix = data_use_MT, col = col_fun, name = "Expression", 
               #row_title = "Signature gene", row_title_side = "right", row_title_rot = 0, row_title_gp = gpar(fontsize = 12), 
               column_title = "Cell", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
               cluster_rows = F, #row_split = des_mat$class, row_gap = unit(0.3, "mm"), border = T, 
               cluster_columns = F, 
               show_row_names = F, row_names_gp = gpar(fontsize = 3), 
               show_column_names = F, column_names_gp = gpar(fontsize = 12), column_names_rot = 45, 
               top_annotation = columnAnnotation(Tissue = cell_metatable_sel$tissue, Group = cell_metatable_sel$group_new, 
                                                 col = list(Tissue = ct_color[as.character(cell_metatable_sel$tissue)], 
                                                            Group = cg_color[as.character(cell_metatable_sel$group_new)]), 
                                                 simple_anno_size = unit(0.3, "cm"), 
                                                 gap = unit(0.08, "cm"), 
                                                 annotation_name_gp = gpar(fontsize = 12)), 
               right_annotation = rowAnnotation(foo = anno_mark(at = match(gene_hlt, rownames(data_use_MT)), labels = gene_hlt, side = "right", labels_gp = gpar(fontsize = 10), link_width = unit(3, "mm"))),
               heatmap_legend_param = list(direction = "vertical", at = c(-2.5, 0, 2.5), legend_height = unit(2.5, "cm"))
)
draw(ht1, row_title = "Signature gene", row_title_gp = gpar(fontsize = 14), 
     merge_legend = F, heatmap_legend_side = "right")
ht_opt(RESET = T)

dev.off()
###

rm(expr_data, expr_ori)
#save.image(file = paste0(OUT, "/Seurat_step2.RData"))
#load(paste0(OUT, "/Seurat_step2.RData"))

save.image(file = paste0(OUT, "/clustering.RData"))
