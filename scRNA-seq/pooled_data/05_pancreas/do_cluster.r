# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data/05_pancreas/")

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

dims_use <- 1:34  # 1:25
resol <- 0.6
expr <- FindClusters(object = expr_ori, reduction.type = "pca", dims.use = dims_use, 
                     resolution = resol, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")
#PrintFindClustersParams(object = expr)
expr@meta.data$cluster <- expr@meta.data[, grep("res.", colnames(expr@meta.data), fixed = T)]
table(expr@meta.data$cluster)

# Run Non-linear dimensional reduction (tSNE)
expr <- RunTSNE(object = expr, dims.use = dims_use, nthreads = 20, do.fast = T)
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T, 
         plot.title = paste0("Dimension: 1:", max(dims_use), " Resolution: ", resol))

### manually curation
xxx <- as.data.frame(expr@dr$tsne@cell.embeddings)

# find another cluster
new_id <- 14
expr@ident <- factor(expr@ident, levels = c(levels(expr@ident), new_id))
xxx_cells <- rownames(subset(xxx, tSNE_1 > -38 & tSNE_1 < -35 & tSNE_2 > -12 & tSNE_2 < -8))
expr@ident[xxx_cells] <- new_id
expr@meta.data[xxx_cells, "cluster"] <- new_id
# find another cluster
new_id <- 15
expr@ident <- factor(expr@ident, levels = c(levels(expr@ident), new_id))
xxx_cells <- intersect(names(expr@ident)[expr@ident == 4], rownames(subset(xxx, tSNE_1 < 17.5)))
expr@ident[xxx_cells] <- new_id
expr@meta.data[xxx_cells, "cluster"] <- new_id
# find another cluster
new_id <- 16
expr@ident <- factor(expr@ident, levels = c(levels(expr@ident), new_id))
xxx_cells <- rownames(subset(xxx, tSNE_1 > -19 & tSNE_1 < -14 & tSNE_2 > -10 & tSNE_2 < -2))
expr@ident[xxx_cells] <- new_id
expr@meta.data[xxx_cells, "cluster"] <- new_id
# find another cluster
new_id <- 17
expr@ident <- factor(expr@ident, levels = c(levels(expr@ident), new_id))
xxx_cells <- rownames(subset(xxx, tSNE_1 > 0 & tSNE_1 < 6 & tSNE_2 > 37.25 & tSNE_2 < 44))
expr@ident[xxx_cells] <- new_id
expr@meta.data[xxx_cells, "cluster"] <- new_id
# find another cluster
new_id <- 18
expr@ident <- factor(expr@ident, levels = c(levels(expr@ident), new_id))
xxx_cells <- rownames(subset(xxx, tSNE_1 > -6 & tSNE_1 < -2 & tSNE_2 > -34 & tSNE_2 < -30))
expr@ident[xxx_cells] <- new_id
expr@meta.data[xxx_cells, "cluster"] <- new_id

rm(xxx, xxx_cells)
###

# check over-clustering
library("foreach")
registerDoSEQ() # avoid connection issue
expr_xxx <- ValidateClusters(expr, pc.use = dims_use, top.genes = 30, min.connectivity = 0.01, acc.cutoff = 0.85, verbose = TRUE)
rm(expr_xxx)

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
cairo_pdf(paste0(OUT, "/Seurat_batchDebug.pdf"), width = 6, height = 6, onefile = T)
do_batchDebug()
dev.off()
###

### check marker (expressed ratio)
marker_genes <- c("PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD14", "CD19")
do_checkMarker(x = expr_data, marker = marker_genes, cell = "")
###

### find markers for every cluster
expr.markers_pn <- do_findMarker(expr, only.pos = F, ncpu = 12)
expr.markers <- subset(expr.markers_pn, avg_logFC > 0)
table(expr.markers$cluster)
expr.markers_ftd <- expr.markers[expr.markers$power>=0.4 & expr.markers$avg_logFC>=log(2), ]
table(expr.markers_ftd$cluster)
top1 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
top3 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(3, avg_logFC)
top5 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10 <- expr.markers_ftd %>% group_by(cluster) %>% top_n(10, avg_logFC)
###

### cell type specific marker genes
expr.markers_spec <- do_findSpecMarker(expr, ncpu = 12)
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
cluster_distr <- acast(data = cellMetaData[, c("bigBatch", "cluster")], formula = bigBatch ~ cluster, fun.aggregate = length, value.var = "cluster")
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
DoHeatmap(object = expr, genes.use = top3$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = expr, genes.use = c(top3$gene, "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA"), slim.col.label = TRUE, remove.key = T)
#DoHeatmap(object = expr, genes.use = c(expr.markers_ftd$gene, "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "IL2RA"), slim.col.label = TRUE, remove.key = T)

top1_rmdup <- as.data.frame(top5)[! duplicated(as.data.frame(top5)$gene), ] %>% group_by(cluster) %>% top_n(1, avg_logFC)
p <- DotPlot_new(object = expr, genes.plot = top1_rmdup$gene, plot.legend = TRUE, x.lab.rot = T, rev.x = T, rev.y = T, do.plot = F, do.return = T)
p + theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))

# genes from literature
#mg_literature <- read.table(paste0(OUT, "/mg_literature.txt"), header = T, sep = "\t", stringsAsFactors = F)

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
enriched_LS <- do_GOenrich(expr.markers_ftd, ncpu = 12)
topEnrich <- do.call("rbind", lapply(enriched_LS, function(x) { y <- x[1:20, c(9,1:8)] }))
View(topEnrich)
###

current.cluster.ids <- as.numeric(levels(expr@ident))
new.cluster.ids <- c("Fibro-PAMR1+SOX6+", "Fibro-PAMR1+SOX6-", "Endo-FLT1", "B-IGHM", "Fibro-IGFBP3", "T", "Fibro-ITGA8", "Fibro-THBS2", "CACNA1A", "B-MME", 
                     "Erythrocyte", "NKT", "Fibro-Prol", "Epi-CADPS", "Endo-NTS", "SM", "Epi-SPINK1", "B-Prol", "Glial")

id_DF <- data.frame(old = as.character(current.cluster.ids), new = new.cluster.ids, stringsAsFactors = F)
id_DF
expr_assigned <- expr
expr_assigned@ident <- plyr::mapvalues(x = expr@ident, from = current.cluster.ids, to = new.cluster.ids)
# change sort
cell_type_sorted <- c("Epi-CADPS", "Epi-SPINK1", "Endo-FLT1", "Endo-NTS", 
                      "SM", "Fibro-PAMR1+SOX6+", "Fibro-PAMR1+SOX6-", "Fibro-THBS2", "Fibro-IGFBP3", "Fibro-ITGA8", "Fibro-Prol", 
                      "Glial", "B-IGHM", "B-MME", "B-Prol", "NKT", "T", "CACNA1A", "Erythrocyte")
cell_type_sorted
expr_assigned@ident <- factor(expr_assigned@ident, levels = cell_type_sorted)

TSNEPlot(object = expr_assigned, do.label = TRUE, pt.size = 1, no.legend = T)

# combine meta with identity
cellMetaDatax <- merge(cellMetaData, data.frame(ident = expr_assigned@ident, stringsAsFactors = F), by.x = 1, by.y = 0, sort = F)
#cellMetaDatax$ident <- factor(cellMetaDatax$ident, levels = sort(levels(cellMetaDatax$ident)))
colnames(cellMetaDatax)[1] <- "cell"
#cellMetaDatax <- cellMetaDatax[! is.na(cellMetaDatax$ident), ]
expr.markers_ftd_labelled <- merge(x = expr.markers_ftd, y = unique(cellMetaDatax[, c("cluster", "ident")]), by = "cluster", sort = F)
expr.markers_spec_LS <- split(expr.markers_spec$gene, expr.markers_spec$cluster)
expr.markers_ftd_labelled$spec <- apply(expr.markers_ftd_labelled[, c("cluster", "gene")], 1, function(x) { 
  x[2] %in% expr.markers_spec_LS[[x[1]]]
})
#expr.markers_ftd_labelled <- expr.markers_ftd_labelled[! is.na(expr.markers_ftd_labelled$ident), ]

# stat by ident
# general correlation
all(colnames(expr@raw.data) == cellMetaDatax$cell)
tmp_LS <- split(as.data.frame(t(expr@raw.data)), cellMetaDatax$ident)

# average cor
cl <- makeCluster(min(length(levels(expr@ident)), 12))
clusterExport(cl = cl, varlist = "tmp_LS")
cor_avg <- parSapply(cl, tmp_LS, function(x) {
  x1 <- t(x)
  y <- sapply(tmp_LS, function(k) {
    x2 <- t(k)
    ky <- mean(cor(x1, x2))
  })
})
stopCluster(cl); rm(cl)

# average profile
tmp <- sapply(tmp_LS, function(x) { y <- colSums(x); y <- y / sum(y); return(y) })

tmp_cor <- cor_avg
tmp_cor[lower.tri(tmp_cor)] <- NA
#corrplot(tmp_cor, type="upper", method = "color", col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
tmp_cor_melted <- na.omit(melt(tmp_cor))
tmp_cor_melted$Var2 <- factor(tmp_cor_melted$Var2, levels = rev(levels(tmp_cor_melted$Var2)))
p1 <- ggplot(tmp_cor_melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradientn(name = "Correlation", colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
  #geom_text(data = subset(tmp_cor_melted, Var1 == Var2), aes(label = Var1), angle = 45, hjust = 0, nudge_x = -1.5, nudge_y = -1.5) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.88))

tmp_cor <- cor(tmp)
tmp_cor[lower.tri(tmp_cor)] <- NA
#corrplot(tmp_cor, type="upper", method = "color", col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
tmp_cor_melted <- na.omit(melt(tmp_cor))
tmp_cor_melted$Var2 <- factor(tmp_cor_melted$Var2, levels = rev(levels(tmp_cor_melted$Var2)))

p2 <- ggplot(tmp_cor_melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradientn(name = "Correlation", colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
  #geom_text(data = subset(tmp_cor_melted, Var1 == Var2), aes(label = Var1), angle = 45, hjust = 0, nudge_x = -1.5, nudge_y = -1.5) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.88))

pdf(file = paste0(OUT, "/", samplingPos, "/Seurat_cellType_heatmap.pdf"), width = 8, height = 8)
print(p1)
print(p2)
dev.off()

# plot GO enrich
enrich_res <- do.call("rbind", enriched_LS)
enrich_res$Term <- Hmisc::capitalize(gsub(" \\(.*\\)", "", enrich_res$Term))
enrich_res$Adjusted.P.value <- enrich_res$q_value
#enrich_res <- enrich_res[order(enrich_res$Adjusted.P.value), ]
enrich_res <- subset(enrich_res, Adjusted.P.value < 0.05)
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res, enrich_res$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, id_DF, by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# # add dummy
nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new), 
                                                 Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
enrich_res_DF$new <- factor(enrich_res_DF$new, levels = cell_type_sorted)
enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", samplingPos, "/enrichment_all.pdf"), width = 8, height = 8)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) + 
  geom_bar(stat = "identity", width = 0.9, show.legend = F) + 
  facet_grid(new ~ ., scales ="free", drop = F) + 
  coord_flip() + 
  ylab(expression(paste(-log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.01)) + 
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") + 
  theme(panel.spacing.y = unit(0.2, "lines"))
gp
dev.off()
#

pdf(paste0(OUT, "/Seurat_tSNE_assigned.pdf"), width = 6, height = 6, useDingbats = F) # 6*6

tmp <- expr_assigned
#tmp@ident <- factor(tmp@ident, levels = sort(levels(tmp@ident)))
ident_labels <- paste(1:length(levels(tmp@ident)), levels(tmp@ident), sep = ": ")
p <- TSNEPlot(object = tmp, do.label = TRUE, pt.size = 1, label.size = 3, no.legend = F, do.return = T) +
  scale_color_manual(labels = ident_labels, values = scales::hue_pal()(length(ident_labels))) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  #guides(color = guide_legend(ncol = 1)) + 
  coord_cartesian(clip = "off") + 
  xlab("tSNE-1") + ylab("tSNE-2")
###
p$layers[[3]]$data[4, 2] <- -33
###
p

# color by seq ID
cellSeqID <- read.table("01-cleandata/merged/cleanFqStat.txt", header = F, sep = "\t", stringsAsFactors = F)[, c(5, 13)]
dim(cellSeqID)
colnames(cellSeqID) <- c("cell", "seqid")
tmp@meta.data$seqid <- cellSeqID[match(rownames(tmp@meta.data), cellSeqID$cell), "seqid"]

DimPlot(object = tmp, reduction.use = "tsne", group.by = "seqid", do.label = F, pt.size = 1, no.legend = F, do.return = T) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank(), legend.background = element_blank()) + 
  theme(legend.position = c(0.135, 0.10)) + 
  xlab("tSNE-1") + ylab("tSNE-2")

# color by bigBatch
if(length(unique(expr@meta.data$bigBatch)) > 1) {
  ### change levels
  tmp@meta.data$bigBatch <- factor(tmp@meta.data$bigBatch, levels = c("YX_C"))
  levels(tmp@meta.data$bigBatch)
  levels(tmp@meta.data$bigBatch) <- c("XXX")
  levels(tmp@meta.data$bigBatch)
  ###
  DimPlot(object = tmp, reduction.use = "tsne", group.by = "bigBatch", do.label = F, pt.size = 1, no.legend = F, do.return = T) + 
    theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank(), legend.background = element_blank()) + 
    theme(legend.position = c(0.125, 0.10)) + 
    xlab("tSNE-1") + ylab("tSNE-2")
}

# hierarchy marker
FeaturePlot_new(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "HBG1"), new.title = c("EPCAM (CD326)", "VIM", "PTPRC (CD45)", "HBG1"), 
                cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 1, no.axes = T)

# dot plot
genes_list <- c("EPCAM", "CADPS", "SPINK1", "PECAM1", "FLT1", "NTS", 
                "ACTA2", "COL1A1", "PAMR1", "SOX6", "THBS2", "IGFBP3", "ITGA8", "TYMS", 
                "PLP1", "PTPRC", "CD79A", "IGHM", "MME", "TOP2A", "FCER1G", "CD3D", 
                "CACNA1A", "HBG1")
p <- DotPlot_new(object = expr_assigned, genes.plot = genes_list, 
                 plot.legend = TRUE, x.lab.rot = T, rev.x = T, rev.y = T, do.plot = F, do.return = T) + 
  theme(legend.position="bottom", legend.justification = "center") + 
  theme(legend.box.margin = margin(r = 50)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5, linetype = 1))
p$layers <- c(geom_rect(xmin = 1-0.5, xmax = 1+0.5, ymin = 18-0.5, ymax = 19+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 4-0.5, xmax = 4+0.5, ymin = 16-0.5, ymax = 17+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 8-0.5, xmax = 8+0.5, ymin = 9-0.5, ymax = 14+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 16-0.5, xmax = 16+0.5, ymin = 3-0.5, ymax = 7+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 17-0.5, xmax = 17+0.5, ymin = 5-0.5, ymax = 7+0.5, fill = "skyblue"), p$layers)
p + guides(fill = guide_colorbar(title = "Z-score", order = 1), size = guide_legend(title = "Fraction (%)", label.vjust = 1.95, label.position = "bottom", override.aes = list(colour = "black")))

# proliferative score (cell cycle)
tmp <- merge(expr_assigned@meta.data, id_DF, by.x = "cluster", by.y = "old", sort = F)
tmp$new <- factor(tmp$new, levels = cell_type_sorted)
ggplot(tmp, aes(x = new, y = S.Score + G2M.Score, fill = new)) + geom_boxplot(show.legend = F) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Proliferative score")

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
  #theme_bw() + 
  facet_grid(. ~ gene, scales = "free_x", switch = "x") + coord_flip() + 
  scale_fill_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  scale_color_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  theme(panel.grid = element_blank()) + 
  theme(strip.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = -15, r = 7.5, b = 15)), strip.background = element_blank()) + 
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
#tmp <- subset(tmp, ! is.na(ident) & spec)
tmp <- tmp[! duplicated(tmp$gene), ]
tmp <- do.call("rbind", split(tmp, tmp$ident))

# violin plot
tmp_data <- data.frame(FetchData(object = expr_assigned, vars.all = unique(tmp$gene)), check.names = F)
tmp_data$cell <- rownames(tmp_data)
tmp_data <- merge(tmp_data, cellMetaDatax[, c("cell", "ident")], by = "cell", sort = F)
tmp_data_melted <- melt(tmp_data, id.vars = c("cell", "ident"), variable.name = "gene")
tmp_data_melted$ident <- factor(tmp_data_melted$ident, levels = rev(levels(tmp_data_melted$ident)))
tmp_data_melted$gene <- factor(tmp_data_melted$gene, levels = unique(tmp$gene))
ggplot(tmp_data_melted, aes(x = ident, y = value, fill = ident, color = ident)) + geom_violin(show.legend = F, scale = "width") + 
  #theme_bw() + 
  facet_grid(. ~ gene, scales = "free_x", switch = "x") + coord_flip() + 
  scale_fill_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  scale_color_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  theme(panel.grid = element_blank()) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5)), strip.background = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linetype = 1, size = 0.75))

# fraction
# tmp <- table(cellMetaDatax[, c("batch", "ident")], useNA = "al")
# tmp_DF <- melt(tmp / rowSums(tmp))
# tmp_DF <- tmp_DF[ (! is.na(tmp_DF$batch)) & (! is.na(tmp_DF$ident) ), ]
# ggplot(tmp_DF, aes(x = ident, y = value * 100, fill = ident)) + 
#   stat_summary(fun.y=mean,position=position_dodge(width=0.95), geom="bar", show.legend = F) + 
#   stat_summary(fun.data=mean_cl_normal, width = 0.6, geom = "errorbar", show.legend = F) + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Cell type") + ylab("Cell type composition (%)") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp_DF$value) * 100 * 1.05))

# fraction (pooled)
tmp <- table(cellMetaDatax[, c("ident")], useNA = "al")
tmp_DF <- melt(tmp / sum(tmp))
colnames(tmp_DF)[1] <- "ident"
tmp_DF <- tmp_DF[ ! is.na(tmp_DF$ident), ]
ggplot(tmp_DF, aes(x = ident, y = value, fill = ident)) + 
  geom_bar(stat = "identity", show.legend = F) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Cell type") + ylab("Cell type composition") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp_DF$value) * 1.05))

# fraction (pooled - pie)
tmp_x <- rep(1.6, nrow(tmp_DF))
tmp_x[8:11] <- seq(1.7, 2, length.out = 4)
ggplot(tmp_DF, aes(x = factor(1), y = value, fill = ident)) + 
  geom_bar(stat = "identity", width = 1, show.legend = F, color = "white") + 
  geom_path(data = data.frame(x = c(1.55, 1.55), y = 1 - c(cumsum(tmp_DF$value)[grep("^Fibro-", tmp_DF$ident)[1] - 1], cumsum(tmp_DF$value)[rev(grep("^Fibro-", tmp_DF$ident))[1]])),
            aes(x = x, y = y), color = "#00BF74", size = 1.5, inherit.aes = F) + 
  geom_path(data = data.frame(x = c(1.55, 1.55), y = 1 - c(cumsum(tmp_DF$value)[12], cumsum(tmp_DF$value)[17])),
            aes(x = x, y = y), color = "mediumpurple1", size = 1.5, inherit.aes = F) + 
  coord_polar(theta = "y", direction = -1) + 
  #geom_text(aes(label = ident, color = ident, x = tmp_x, y = 1 - (cumsum(tmp_DF$value) - tmp_DF$value/2), angle = -(cumsum(tmp_DF$value) - tmp_DF$value/2) * 360), show.legend = F) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
rm(tmp_x)

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

# fraction in different pos (bar with connected lines)
tmp <- table(cellMetaDatax[, c("bigBatch", "ident")])
tmp_DF <- melt(tmp / rowSums(tmp))
# sort by real position
tmp_DF$bigBatch <- factor(tmp_DF$bigBatch, levels = c("YX_C"))

### change levels
levels(tmp_DF$bigBatch)
levels(tmp_DF$bigBatch) <- c("mixture")
levels(tmp_DF$bigBatch)

tmp_DF_sub <- subset(tmp_DF, bigBatch == rev(levels(tmp_DF$bigBatch))[1])

ggplot(tmp_DF, aes(x = bigBatch, y = value, fill = ident)) + geom_bar(stat = "identity", width = 0.8, color = "white") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0, 0)) + 
  #geom_segment(mapping = aes(x = 1+0.8/2, xend = 2-0.8/2, y = 1, yend = 1), size = 0.3) + 
  #geom_segment(data = dcast(tmp_DF, ident ~ bigBatch), mapping = aes(x = 1+0.8/2, xend = 2-0.8/2, y = (1-cumsum(Central)), yend = (1-cumsum(Surrounding))), size = 0.3) + 
  geom_path(data = data.frame(x = c(1.5, 1.5), y = 1 - c(cumsum(tmp_DF_sub$value)[grep("^Fibro-", tmp_DF_sub$ident)[1] - 1], cumsum(tmp_DF_sub$value)[rev(grep("^Fibro-", tmp_DF_sub$ident))[1]])), 
            aes(x = x, y = y), color = "#00BF74", size = 1.5, inherit.aes = F) + 
  geom_path(data = data.frame(x = c(1.5, 1.5), y = 1 - c(cumsum(tmp_DF_sub$value)[12], cumsum(tmp_DF_sub$value)[17])),
            aes(x = x, y = y), color = "mediumpurple1", size = 1.5, inherit.aes = F) + 
  guides(fill = guide_legend(ncol = 1)) + 
  ylab("Cell type composition")

dev.off()

###
pdf(paste0(OUT, "/Seurat_tSNE_assigned_big.pdf"), width = 8, height = 6, useDingbats = F)

# dot plot
genes_list <- c("EPCAM", "CADPS", "SPINK1", "PECAM1", "FLT1", "NTS", 
                "ACTA2", "COL1A1", "PAMR1", "SOX6", "THBS2", "IGFBP3", "ITGA8", "TYMS", 
                "PLP1", "PTPRC", "CD79A", "IGHM", "MME", "TOP2A", "FCER1G", "CD3D", 
                "CACNA1A", "HBG1")
p <- DotPlot_new(object = expr_assigned, genes.plot = genes_list, 
                 plot.legend = TRUE, x.lab.rot = T, rev.x = T, rev.y = T, do.plot = F, do.return = T) + 
  theme(legend.position="bottom", legend.justification = "center") + 
  theme(legend.box.margin = margin(r = 50)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5, linetype = 1))
p$layers <- c(geom_rect(xmin = 1-0.5, xmax = 1+0.5, ymin = 18-0.5, ymax = 19+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 4-0.5, xmax = 4+0.5, ymin = 16-0.5, ymax = 17+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 8-0.5, xmax = 8+0.5, ymin = 9-0.5, ymax = 14+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 16-0.5, xmax = 16+0.5, ymin = 3-0.5, ymax = 7+0.5, fill = "skyblue"), p$layers)
p$layers <- c(geom_rect(xmin = 17-0.5, xmax = 17+0.5, ymin = 5-0.5, ymax = 7+0.5, fill = "skyblue"), p$layers)
p + guides(fill = guide_colorbar(title = "Z-score", order = 1), size = guide_legend(title = "Fraction (%)", label.vjust = 1.95, label.position = "bottom", override.aes = list(colour = "black")))

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
  #theme_bw() + 
  facet_grid(. ~ gene, scales = "free_x", switch = "x") + coord_flip() + 
  scale_fill_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  scale_color_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  theme(panel.grid = element_blank()) + 
  theme(strip.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = -15, r = 7.5, b = 15)), strip.background = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linetype = 1, size = 0.75))
pg <- ggplotGrob(p)
for(i in which(grepl("strip-b", pg$layout$name))){ pg$grobs[[i]]$layout$clip <- "off" }
grid.newpage(); grid::grid.draw(pg)

# collagen
tmp <- expr.markers_ftd_labelled[grep("^COL[0-9]", expr.markers_ftd_labelled$gene, value = F), ]
tmp <- tmp[order(as.character(tmp$ident)), ]
#tmp <- subset(tmp, ! is.na(ident) & spec)
tmp <- tmp[! duplicated(tmp$gene), ]
tmp <- do.call("rbind", split(tmp, tmp$ident))

# violin plot
tmp_data <- data.frame(FetchData(object = expr_assigned, vars.all = unique(tmp$gene)), check.names = F)
tmp_data$cell <- rownames(tmp_data)
tmp_data <- merge(tmp_data, cellMetaDatax[, c("cell", "ident")], by = "cell", sort = F)
tmp_data_melted <- melt(tmp_data, id.vars = c("cell", "ident"), variable.name = "gene")
tmp_data_melted$ident <- factor(tmp_data_melted$ident, levels = rev(levels(tmp_data_melted$ident)))
tmp_data_melted$gene <- factor(tmp_data_melted$gene, levels = unique(tmp$gene))
ggplot(tmp_data_melted, aes(x = ident, y = value, fill = ident, color = ident)) + geom_violin(show.legend = F, scale = "width") + 
  #theme_bw() + 
  facet_grid(. ~ gene, scales = "free_x", switch = "x") + coord_flip() + 
  scale_fill_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  scale_color_manual(values = rev(scales::hue_pal()(length(levels(tmp_data_melted$ident))))) + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) + 
  theme(panel.grid = element_blank()) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5)), strip.background = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linetype = 1, size = 0.75))

dev.off()
###

pdf(paste0(OUT, "/Seurat_tSNE_assigned_heatmap.pdf"), width = 12, height = 12, useDingbats = F)
# marker
tmp <- expr.markers_ftd_labelled %>% group_by(cluster) %>% top_n(6, avg_logFC)
tmp <- tmp[order(as.character(tmp$ident)), ]
tmp <- subset(tmp, ! is.na(ident))
tmp <- tmp[! duplicated(tmp$gene), ]
tmp <- do.call("rbind", split(tmp, tmp$ident))
DoHeatmap_new(object = expr_assigned, genes.use = tmp$gene, genes.group = tmp$ident, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, do.colBar = T, strip.text.x.top = 15, 
              panel.spacing.y = 0, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA), max.cells.per.ident = 80, random.seed = 1))

# cell type specific marker
tmp <- subset(expr.markers_ftd_labelled, spec) %>% group_by(cluster) %>% top_n(6, avg_logFC)
tmp <- tmp[order(as.character(tmp$ident)), ]
tmp$ident <- factor(tmp$ident)
tmp <- subset(tmp, ! is.na(ident))
tmp <- do.call("rbind", split(tmp, tmp$ident))
DoHeatmap_new(object = expr_assigned, genes.use = tmp$gene, genes.group = tmp$ident, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, do.colBar = T, strip.text.x.top = 15, 
              panel.spacing.y = 0, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA), max.cells.per.ident = 80, random.seed = 1))

# literature
# tmp <- do.call("rbind", lapply(split(mg_literature, mg_literature$Cluster), function(x) head(x, n = 5)))
# DoHeatmap_new(object = expr_assigned, genes.use = tmp$gene, genes.group = factor(tmp$Cluster), slim.col.label = TRUE, remove.key = F, rotate.key = T, 
#               group.cex = 11, cex.row = 7, do.colBar = T, 
#               group.order = sort(levels(expr_assigned@ident)), panel.spacing.y = 0, 
#               cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA))) + 
#   theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 8, r = 0, b = 0, l = 0)), strip.text.y = element_blank()) + 
#   theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 0, unit = "cm"))

# cell cycle
cyc_DF <- data.frame(gene = c(s.genes, g2m.genes), stage = rep(c("S", "G2/M"), c(length(s.genes), length(g2m.genes))), stringsAsFactors = F)
cyc_DF_sub <- subset(cyc_DF, gene %in% rownames(expr_assigned@data))
tmp_data <- data.frame(FetchData(object = expr_assigned, vars.all = cyc_DF_sub$gene), check.names = F)
tmp_DF <- data.frame(expr = rowMeans(tmp_data), stringsAsFactors = F)
tmp_DF <- merge(tmp_DF, cellMetaDatax, by.x = 0, by.y = "cell", sort = F)
colnames(tmp_DF)[1] <- "cell"
tmp_DF_sorted <- tmp_DF[order(tmp_DF$ident, - tmp_DF$expr), ]
cols1 <- scales::hue_pal()(length(levels(tmp_DF_sorted$ident)))
cols2 <- c("chartreuse", "lightblue")

# p <- DoHeatmap_new(object = expr_assigned, genes.use = cyc_DF$gene, genes.group = cyc_DF$stage, slim.col.label = TRUE, remove.key = F, rotate.key = T,
#                    group.cex = 11, cex.row = 7, switch_y = T,
#                    col.low = "blue", col.mid = "white", col.high = "red",
#                    group.order = sort(levels(expr_assigned@ident)), panel.spacing.y = 0.003,
#                    cells.use = tmp_DF_sorted$cell) +
#   theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 1, l = 0)), strip.text.y = element_text(size = 11, margin = margin(t = 0, r = 5, b = 0, l = 5))) + 
#   theme(strip.background.x = element_rect(color = NA), strip.background.y = element_rect(color = NA)) + 
#   theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = 0, unit = "cm")) +
#   theme(axis.text.y = element_text(size = 4))
# 
# g <- ggplot_gtable(ggplot_build(p))
# strip_both <- which(grepl('strip-', g$layout$name))
# fills <- c(cols1, cols2)
# k <- 1
# for (i in strip_both) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.newpage()
# grid.draw(g)

DoHeatmap_new(object = expr_assigned, genes.use = cyc_DF$gene, genes.group = factor(cyc_DF$stage, levels = c("S", "G2/M")), slim.col.label = TRUE, remove.key = F, rotate.key = T, 
              group.cex = 11, cex.row = 7, do.colBar = T, switch_y = T, strip.text.x.top = 15, 
              col.low = "blue", col.mid = "white", col.high = "red",
              #group.order = sort(levels(expr_assigned@ident)), 
              panel.spacing.y = 0.01, 
              cells.use = WhichCells(expr_assigned, ident = setdiff(levels(expr_assigned@ident), NA), max.cells.per.ident = 80, random.seed = 1), 
              strip.text.y.display = T, strip.text.y = 12)
dev.off()

# check marker after assign identity
pdf(paste0(OUT, "/Seurat_tSNE_spMarker.pdf"), width = 6, height = 6, useDingbats = F)
DoHeatmap(object = expr_assigned, genes.use = c("EPCAM", "CD3D", "CD4", "CD8A", "HBA1"), 
          slim.col.label = TRUE, remove.key = TRUE, cex.row = 8) + 
  theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5, r = 0, b = 0, l = 0)))

mg_DF <- data.frame(ident = expr_assigned@ident, expressed = expr@raw.data["EPCAM", names(expr@ident), drop = T] > 0, stringsAsFactors = F)
#mg_DF$ident <- factor(mg_DF$ident, levels = sort(levels(mg_DF$ident)))
ggplot(mg_DF, aes(x = ident, fill = expressed)) + geom_bar() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(table(mg_DF$ident)*1.05))) + 
  scale_fill_discrete(name = "EPCAM") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + 
  theme(legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + ylab("Cell number")

# FeaturePlot
for(i in c("EPCAM", "VIM", "PECAM1", "PROM1", "COL1A1", "IRF1", "FOSB", "CACNA1A", "TYMS", top1$gene)) {
  p <- FeaturePlot_new(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, do.plot = F, do.return = T)[[1]]
  p <- p + theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) +
    xlab("tSNE-1") + ylab("tSNE-2")
  print(p)
}
dev.off()

# 4. save and write meta table ----
#save(expr, file = paste0(OUT, "/Seurat_expr.Robj"))
write.table(x = cellMetaDatax, file = paste0(OUT, "/Seurat_metaData.txt"), row.names = F, col.names = T, quote = F,sep = "\t")
write.table(x = expr.markers_pn, file = paste0(OUT, "/Seurat_markerGenes_pn.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = expr.markers_ftd_labelled, file = paste0(OUT, "/Seurat_markerGenes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
data.table::fwrite(x = as.data.frame(expr@scale.data), file = paste0(OUT, "/UMIcount_scaled.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)

save.image(file = paste0(OUT, "/clustering.RData"))
