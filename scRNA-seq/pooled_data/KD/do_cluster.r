# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data/KD/")

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

# Load module genes
md_genes <- read.table("../../pooled_data/All/03-expression/merged/geneModule/geneModule_genes.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(md_genes) <- c("mdid", "gene")

# Load the dataset
expr_data <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_cellFiltered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
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
expr <- JackStrawPlot(object = expr, PCs = 1:50)
PCElbowPlot(object = expr, num.pc = 50)

dim_limit <- min(which(expr@dr$pca@jackstraw@overall.p.values[, 2] > 0.001)) - 1
cat("Suggested max. dimension:", dim_limit, "\n")

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

dims_use <- 1:dim_limit  # 1:25
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

DimPlot(object = expr, reduction.use = "tsne", group.by = "batch", no.legend = F, do.return = TRUE, 
        vector.friendly = TRUE, pt.size = 3) + ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

### manually curation
xxx <- as.data.frame(expr@dr$tsne@cell.embeddings)

# cell type in 20w
ctMeta_20w <- read.table(file = "../../pooled_data/All/cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
ts_id <- "small intestine"
ctMeta_20w <- subset(ctMeta_20w, tissue == ts_id)
cellType_20w <- ctMeta_20w$ident

# use the ident from CellBlast
cell2ident <- read.table("cell2ident.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cell2ident)
colnames(cell2ident) <- c("cell", "ident")
cell2ident <- cell2ident[match(names(expr@ident), cell2ident$cell), ]
as.data.frame(table(cell2ident$ident))
cell2ident$ori.ident <- cell2ident$ident
cell2ident$ident[cell2ident$ident %in% c("ambiguous", "rejected")] <- "Unknown"

expr_tmp <- expr
expr_tmp@ident <- factor(cell2ident$ident, levels = c(cellType_20w, setdiff(unique(cell2ident$ident), cellType_20w)))
names(expr_tmp@ident) <- cell2ident$cell
expr_tmp@meta.data$cluster <- cell2ident$ident
gp <- TSNEPlot(object = expr_tmp, pt.size = 2, do.label = T, no.legend = T, colors.use = c(ctMeta_20w[ctMeta_20w$ident %in% unique(expr_tmp@ident), "color"], "#B3B3B3"), do.return = T)
#ggsave(filename = paste0(OUT, "/Seurat_tSNE_before_KNN.pdf"), plot = gp, width = 6, height = 6, useDingbats = F)
rm(expr_tmp)

# rescue ambiguous ident using KNN
cell2ident_known <- subset(cell2ident, ori.ident != "ambiguous" & ori.ident != "rejected")
cell2ident_ambiguous <- subset(cell2ident, ori.ident == "ambiguous")
pca_embeddings <- expr@dr$pca@cell.embeddings[, 1:dim_limit]
knn_res <- class::knn(train = pca_embeddings[cell2ident_known$cell, ], test = pca_embeddings[cell2ident_ambiguous$cell, ], cl = cell2ident_known$ident, k = 10, prob = F)
knn_res <- as.character(knn_res)
cell2ident[match(cell2ident_ambiguous$cell, cell2ident$cell), "ident"] <- knn_res
cell2ident$knn.ident <- cell2ident$ident

expr@ident <- factor(cell2ident$ident, levels = c(cellType_20w, "Unknown"))
names(expr@ident) <- cell2ident$cell
expr@meta.data$cluster <- cell2ident$ident
gp <- TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T, colors.use = c(ctMeta_20w[ctMeta_20w$ident %in% unique(expr@ident), "color"], "#B3B3B3"), do.return = T)
#ggsave(filename = paste0(OUT, "/Seurat_tSNE_after_KNN.pdf"), plot = gp, width = 6, height = 6, useDingbats = F)

# correct ident
cell2ident[cell2ident$cell %in% c("XC_NR2F1-2112_13", "XC_NR2F1-2922_89", "XC_NCFAM-001_96"), "ident"] <- "Unknown"

expr@ident <- factor(cell2ident$ident, levels = c(intersect(cellType_20w, unique(expr@ident)), "Unknown"))
names(expr@ident) <- cell2ident$cell
expr@meta.data$cluster <- cell2ident$ident
gp <- TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T, colors.use = c(ctMeta_20w[ctMeta_20w$ident %in% unique(expr@ident), "color"], "#B3B3B3"), do.return = T)
#ggsave(filename = paste0(OUT, "/Seurat_tSNE_after_correction1.pdf"), plot = gp, width = 6, height = 6, useDingbats = F)

# load signature genes
# snGene <- read.table("../../pooled_data/02_small_intestine/03-expression/merged/cellCluster/Seurat_markerGenes.txt", header = T, sep = "\t", stringsAsFactors = F)
# snGene_top1 <- do.call("rbind", lapply(split(snGene, snGene$ident), function(x) { head(x, 1) }))

# pdf(file = file.path(OUT, "snGene.pdf"), width = 6, height = 6)
# apply(subset(snGene_top1, gene %in% rownames(expr@data), c("ident", "gene")), 1, function(x) {
#   FeaturePlot(object = expr, features.plot = x[2], cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.return = T) #+ ggtitle(x[1])
# })
# dev.off()

# find another cluster
# new_id <- 18
# expr@ident <- factor(expr@ident, levels = c(levels(expr@ident), new_id))
# xxx_cells <- rownames(subset(xxx, tSNE_1 > 14 & tSNE_1 < 15 & tSNE_2 > 35 & tSNE_2 < 37))
# expr@ident[xxx_cells] <- new_id
# expr@meta.data[xxx_cells, "cluster"] <- new_id

# add cond info
expr@meta.data$cond <- gsub("^XC_", "", expr@meta.data$bigBatch)
expr@meta.data$cond <- factor(expr@meta.data$cond, levels = unique(expr@meta.data$cond))
#
levels(expr@meta.data$cond)
levels(expr@meta.data$cond) <- c("Control", "siFOXL1", "siNR2F1")
levels(expr@meta.data$cond)
#
table(expr@meta.data$cond)
expr@meta.data$cond_cluster <- paste(expr@meta.data$cond, expr@meta.data$cluster, sep = ".")

save.image(file = paste0(OUT, "/Seurat_step2.RData"))
#load(paste0(OUT, "/Seurat_step2.RData"))

# 3.1 KD FOXL1 ----

### remove isZero (Control) and nonZero (si)
cell_isZero <- intersect(rownames(expr@meta.data)[expr@meta.data$cond == "Control"], colnames(expr@raw.data)[expr@raw.data["FOXL1", ] <= 2])
expr@meta.data[cell_isZero, "cond_cluster"] <- paste(expr@meta.data[cell_isZero, "cond_cluster"], "isZero", sep = ".")
table(subset(expr@meta.data, grepl("isZero", cond_cluster), "cond", drop = T))
cell_nonZero <- intersect(rownames(expr@meta.data)[expr@meta.data$cond == "siFOXL1"], colnames(expr@raw.data)[expr@raw.data["FOXL1", ] > 0])
expr@meta.data[cell_nonZero, "cond_cluster"] <- paste(expr@meta.data[cell_nonZero, "cond_cluster"], "nonZero", sep = ".")
table(subset(expr@meta.data, grepl("nonZero", cond_cluster), "cond", drop = T))
###
ct_tmp <- apply(data.table::CJ(levels(expr@meta.data$cond), levels(expr@ident), sorted = F), 1, function(x) { paste(x, collapse = ".") })
expr@ident <- factor(expr@meta.data$cond_cluster, levels = c(intersect(ct_tmp, expr@meta.data$cond_cluster), setdiff(unique(expr@meta.data$cond_cluster), ct_tmp)))
names(expr@ident) <- cell2ident$cell
expr@meta.data$cluster_pure <- expr@meta.data$cluster
expr@meta.data$cluster <- expr@ident
#

rm(xxx)
###

# check over-clustering
# library("foreach")
# registerDoSEQ() # avoid connection issue
# expr_xxx <- ValidateClusters(expr, pc.use = dims_use, top.genes = 30, min.connectivity = 0.01, acc.cutoff = 0.85, verbose = TRUE)
# rm(expr_xxx)

# UMAP
expr <- RunUMAP(object = expr, reduction.use = "pca", dims.use = dims_use, min_dist = 1)
#DimPlot(object = expr, reduction.use = "umap", no.legend = F, do.return = TRUE, 
#        vector.friendly = TRUE, pt.size = 3) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

# add info to meta
cellMetaData <- merge(expr@meta.data, expr@dr$tsne@cell.embeddings, by = 0, sort = F)
cellMetaData$batch <- factor(cellMetaData$batch, levels = unique(cellMetaData$batch))
#cellMetaData$cluster <- factor(cellMetaData$cluster, levels = sort(as.numeric(unique(cellMetaData$cluster))))
cellStat <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/filtering_cells.txt"), header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
cellMetaData <- merge(cellMetaData, cellStat[, -c(4, 5)], by.x = 1, by.y = 0, sort = F)
colnames(cellMetaData)[1] <- "cell"

case_TF <- "FOXL1"
ct_1 <- grep("^Control.*COL6A5$", levels(expr@ident), value = T)
ct_2 <- grep("^siFOXL1.*COL6A5$", levels(expr@ident), value = T)
expr.markers_1a <- FindMarkers(object = expr, ident.1 = ct_1, ident.2 = ct_2, test.use = "roc", only.pos = T, min.pct = 0.25)
expr.markers_1a$avg_logFC <- - expr.markers_1a$avg_logFC
expr.markers_1b <- FindMarkers(object = expr, ident.2 = ct_2, ident.1 = ct_1, test.use = "roc", only.pos = T, min.pct = 0.25)

case_mdid <- "MD91"
expr.markers_1s <- rbind(expr.markers_1a, NULL)
expr.markers_1s$gene <- rownames(expr.markers_1s)
expr.markers_1s$isInMd <- expr.markers_1s$gene %in% subset(md_genes, mdid %in% case_mdid, "gene", drop = T)
table(expr.markers_1s$isInMd)
expr.markers_1s_ftd <- subset(expr.markers_1s, power >= 0.4)
table(expr.markers_1s_ftd$isInMd)
expr.markers_1s_ftd <- expr.markers_1s_ftd[order(expr.markers_1s_ftd$avg_logFC, decreasing = T), ]
expr.markers_1s_ftd$rank <- 1:nrow(expr.markers_1s_ftd)

# write (1a + 1b)
expr.markers_1t <- rbind(expr.markers_1a, expr.markers_1b)
expr.markers_1t$gene <- rownames(expr.markers_1t)
expr.markers_1t$isInMd <- expr.markers_1t$gene %in% subset(md_genes, mdid %in% case_mdid, "gene", drop = T)
table(expr.markers_1t$isInMd)
expr.markers_1t_ftd <- subset(expr.markers_1t, power >= 0.4)
table(expr.markers_1t_ftd$isInMd)
expr.markers_1t_ftd <- expr.markers_1t_ftd[order(expr.markers_1t_ftd$avg_logFC, decreasing = F), ]
#expr.markers_1t_ftd$rank <- 1:nrow(expr.markers_1t_ftd)
write.table(x = expr.markers_1t_ftd, file = file.path(OUT, "Seurat_markerGenes_Control_siFOXL1.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#

# significance
expr.markers_1s_tmp <- subset(expr.markers_1s_ftd, gene != case_TF)
sig_pvalue <- phyper(q = sum(expr.markers_1s_tmp$isInMd), m = table(md_genes$mdid)[case_mdid], n = nrow(md_genes) - table(md_genes$mdid)[case_mdid], k = nrow(expr.markers_1s_tmp), lower.tail = F, log.p = FALSE)
sig_pvalue

# plot
pdf(file = paste0(OUT, "/geneExpr_barplot_si", case_TF, ".pdf"), width = 6, height = 4.5, useDingbats = F)
ggplot(expr.markers_1s_ftd, aes(x = rank, y = avg_logFC, fill = isInMd)) + geom_col(show.legend = F) + xlab("Rank") + ylab("Log (siFOXL1 / Control)") + 
  scale_fill_manual(values = c("grey", "tomato")) + theme(aspect.ratio = 0.4) + 
  annotate(geom = "text", x = 0, y = min(expr.markers_1s_ftd$avg_logFC * 0.95), label = paste("P-value:", signif(sig_pvalue, 3)), hjust = 0)
dev.off()

pdf(file = paste0(OUT, "/geneExpr_heatmap_si", case_TF, "_", case_mdid, ".pdf"), width = 6, height = 4.5, useDingbats = F, onefile = F)
###
expr_plot <- expr
ct_1s <- gsub("\\..*", "", levels(expr_plot@ident)[levels(expr_plot@ident) == ct_1])
ct_2s <- gsub("\\..*", "", levels(expr_plot@ident)[levels(expr_plot@ident) == ct_2])
levels(expr_plot@ident)[levels(expr_plot@ident) == ct_1] <- ct_1s
levels(expr_plot@ident)[levels(expr_plot@ident) == ct_2] <- ct_2s
expr_plot@meta.data$cluster <- expr_plot@ident
###
expr.markers_1s_sub <- subset(expr.markers_1s_ftd, isInMd)
gene_tmp <- data.frame(gene = unique(c(case_TF, expr.markers_1s_sub$gene[order(expr.markers_1s_sub$myAUC, decreasing = T)])), stringsAsFactors = F)
gene_tmp$group <- factor(c("TF", rep("gene", nrow(gene_tmp) - 1)), levels = c("TF", "gene"))
DoHeatmap_new(object = expr_plot, genes.use = gene_tmp$gene, genes.group = gene_tmp$group, 
          col.low = "blue", col.mid = "white", col.high = "red", 
          do.colBar = T, colBar.y = 0.8, colBar.col = c("grey", "red"), 
          cells.use = WhichCells(expr_plot, ident = c(ct_1s, ct_2s), max.cells.per.ident = 20, random.seed = 1), 
          slim.col.label = T, remove.key = F, group.label.rot = F, rotate.key = T, panel.spacing.y = 0, 
          group.cex = 12, group.spacing = 0.5, strip.text.x.top = 15, legend.margin.for.colBar = margin(t = -5))
rm(expr_plot)
dev.off()

# GO enrichment
expr.markers_ftd_updw <- expr.markers_1s_sub
expr.markers_ftd_updw <- subset(expr.markers_ftd_updw, gene != case_TF)
expr.markers_ftd_updw$cluster <- ifelse(expr.markers_ftd_updw$avg_logFC > 0, "up", "down")
enriched_LS <- do_GOenrich(expr.markers_ftd_updw, ncpu = 2)

# plot GO enrich
enrich_res <- do.call("rbind", enriched_LS)
enrich_res$Term <- Hmisc::capitalize(gsub(" \\(.*\\)", "", enrich_res$Term))
enrich_res$Adjusted.P.value <- enrich_res$q_value
#enrich_res <- enrich_res[order(enrich_res$Adjusted.P.value), ]
enrich_res <- subset(enrich_res, Adjusted.P.value < 0.05)
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res, enrich_res$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, data.frame(old = c("up", "down"), new = c("up", "down"), stringsAsFactors = F), by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# add dummy
# nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
# enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new),
#                                                  Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
# enrich_res_DF$new <- factor(enrich_res_DF$new, levels = cell_type_sorted)
# enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", samplingPos, "/enrichment_all.pdf"), width = 6, height = 3)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) +
  geom_bar(stat = "identity", width = 0.9, show.legend = F) +
  #facet_grid(new ~ ., scales ="free", drop = F) +
  coord_flip() +
  ylab(expression(paste(-Log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.2)) +
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") +
  theme(panel.spacing.y = unit(0.2, "lines"))
gp
dev.off()
#


# 3.2 KD NR2F1 ----
case_TF <- "NR2F1"
case_mdid <- "MD203"

### remove isZero (Control) and nonZero (si)
cell_isZero <- intersect(rownames(expr@meta.data)[expr@meta.data$cond == "Control"], colnames(expr@raw.data)[expr@raw.data[case_TF, ] <= 2])
expr@meta.data[cell_isZero, "cond_cluster"] <- paste(expr@meta.data[cell_isZero, "cond_cluster"], "isZero", sep = ".")
table(subset(expr@meta.data, grepl("isZero", cond_cluster), "cond", drop = T))
cell_nonZero <- intersect(rownames(expr@meta.data)[expr@meta.data$cond == paste0("si", case_TF)], colnames(expr@raw.data)[expr@raw.data[case_TF, ] > 0])
expr@meta.data[cell_nonZero, "cond_cluster"] <- paste(expr@meta.data[cell_nonZero, "cond_cluster"], "nonZero", sep = ".")
table(subset(expr@meta.data, grepl("nonZero", cond_cluster), "cond", drop = T))
###
ct_tmp <- apply(data.table::CJ(levels(expr@meta.data$cond), levels(expr@ident), sorted = F), 1, function(x) { paste(x, collapse = ".") })
expr@ident <- factor(expr@meta.data$cond_cluster, levels = c(intersect(ct_tmp, expr@meta.data$cond_cluster), setdiff(unique(expr@meta.data$cond_cluster), ct_tmp)))
names(expr@ident) <- cell2ident$cell
expr@meta.data$cluster_pure <- expr@meta.data$cluster
expr@meta.data$cluster <- expr@ident

# add info to meta
cellMetaData <- merge(expr@meta.data, expr@dr$tsne@cell.embeddings, by = 0, sort = F)
cellMetaData$batch <- factor(cellMetaData$batch, levels = unique(cellMetaData$batch))
#cellMetaData$cluster <- factor(cellMetaData$cluster, levels = sort(as.numeric(unique(cellMetaData$cluster))))
cellStat <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/filtering_cells.txt"), header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cellStat)
cellMetaData <- merge(cellMetaData, cellStat[, -c(4, 5)], by.x = 1, by.y = 0, sort = F)
colnames(cellMetaData)[1] <- "cell"

ct_1 <- grep("^Control.*COL6A5$", levels(expr@ident), value = T)
ct_2 <- grep(paste0("^si", case_TF, ".*COL6A5$"), levels(expr@ident), value = T)
expr.markers_1a <- FindMarkers(object = expr, ident.1 = ct_1, ident.2 = ct_2, test.use = "roc", only.pos = T, min.pct = 0.25)
expr.markers_1a$avg_logFC <- - expr.markers_1a$avg_logFC
expr.markers_1b <- FindMarkers(object = expr, ident.2 = ct_2, ident.1 = ct_1, test.use = "roc", only.pos = T, min.pct = 0.25)

expr.markers_1s <- rbind(expr.markers_1a, NULL)
expr.markers_1s$gene <- rownames(expr.markers_1s)
expr.markers_1s$isInMd <- expr.markers_1s$gene %in% subset(md_genes, mdid %in% case_mdid, "gene", drop = T)
table(expr.markers_1s$isInMd)
expr.markers_1s_ftd <- subset(expr.markers_1s, power >= 0.3)
table(expr.markers_1s_ftd$isInMd)
expr.markers_1s_ftd <- expr.markers_1s_ftd[order(expr.markers_1s_ftd$avg_logFC, decreasing = T), ]
expr.markers_1s_ftd$rank <- 1:nrow(expr.markers_1s_ftd)

# write (1a + 1b)
expr.markers_1t <- rbind(expr.markers_1a, expr.markers_1b)
expr.markers_1t$gene <- rownames(expr.markers_1t)
expr.markers_1t$isInMd <- expr.markers_1t$gene %in% subset(md_genes, mdid %in% case_mdid, "gene", drop = T)
table(expr.markers_1t$isInMd)
expr.markers_1t_ftd <- subset(expr.markers_1t, power >= 0.3)
table(expr.markers_1t_ftd$isInMd)
expr.markers_1t_ftd <- expr.markers_1t_ftd[order(expr.markers_1t_ftd$avg_logFC, decreasing = F), ]
#expr.markers_1t_ftd$rank <- 1:nrow(expr.markers_1t_ftd)
write.table(x = expr.markers_1t_ftd, file = file.path(OUT, paste0("Seurat_markerGenes_Control_si", case_TF, ".txt")), row.names = F, col.names = T, quote = F, sep = "\t")
#

# significance
expr.markers_1s_tmp <- subset(expr.markers_1s_ftd, gene != case_TF)
sig_pvalue <- phyper(q = sum(expr.markers_1s_tmp$isInMd), m = table(md_genes$mdid)[case_mdid], n = nrow(md_genes) - table(md_genes$mdid)[case_mdid], k = nrow(expr.markers_1s_tmp), lower.tail = F, log.p = FALSE)
sig_pvalue

# plot
pdf(file = paste0(OUT, "/geneExpr_barplot_si", case_TF, ".pdf"), width = 6, height = 4.5, useDingbats = F)
ggplot(expr.markers_1s_ftd, aes(x = rank, y = avg_logFC, fill = isInMd)) + geom_col(show.legend = F) + xlab("Rank") + ylab(paste0("Log (si", case_TF, " / Control)")) + 
  scale_fill_manual(values = c("grey", "tomato")) + theme(aspect.ratio = 0.4) + 
  annotate(geom = "text", x = 0, y = min(expr.markers_1s_ftd$avg_logFC * 0.95), label = paste("P-value:", signif(sig_pvalue, 3)), hjust = 0)
dev.off()

pdf(file = paste0(OUT, "/geneExpr_heatmap_si", case_TF, "_", case_mdid, ".pdf"), width = 6, height = 4.5, useDingbats = F, onefile = F)
###
expr_plot <- expr
ct_1s <- gsub("\\..*", "", levels(expr_plot@ident)[levels(expr_plot@ident) == ct_1])
ct_2s <- gsub("\\..*", "", levels(expr_plot@ident)[levels(expr_plot@ident) == ct_2])
levels(expr_plot@ident)[levels(expr_plot@ident) == ct_1] <- ct_1s
levels(expr_plot@ident)[levels(expr_plot@ident) == ct_2] <- ct_2s
expr_plot@meta.data$cluster <- expr_plot@ident
###
expr.markers_1s_sub <- subset(expr.markers_1s_ftd, isInMd)
gene_tmp <- data.frame(gene = unique(c(case_TF, expr.markers_1s_sub$gene[order(expr.markers_1s_sub$myAUC, decreasing = T)])), stringsAsFactors = F)
gene_tmp$group <- factor(rep(tail(levels(expr_plot@ident), 1), nrow(gene_tmp)), levels = levels(expr_plot@ident))
DoHeatmap_new(object = expr_plot, genes.use = gene_tmp$gene, genes.group = gene_tmp$group, 
              col.low = "blue", col.mid = "white", col.high = "red", 
              do.colBar = T, colBar.y = 0.82, colBar.col = c("grey", "red"), 
              cells.use = WhichCells(expr_plot, ident = c(ct_1s, ct_2s), max.cells.per.ident = 30, random.seed = 1), 
              slim.col.label = T, remove.key = F, group.label.rot = F, rotate.key = T, 
              group.cex = 12, group.spacing = 0.5, strip.text.x.top = 15, legend.margin.for.colBar = margin(t = -5))
rm(expr_plot)
dev.off()

### 3. Assigning cell type identity to clusters ----
expr_assigned <- expr

# combine meta with identity
cellMetaDatax <- merge(cellMetaData, data.frame(ident = expr_assigned@ident, stringsAsFactors = F), by.x = 1, by.y = 0, sort = F)
#cellMetaDatax$ident <- factor(cellMetaDatax$ident, levels = sort(levels(cellMetaDatax$ident)))
colnames(cellMetaDatax)[1] <- "cell"
#cellMetaDatax <- cellMetaDatax[! is.na(cellMetaDatax$ident), ]

# 4. save and write meta table ----
#save(expr, file = paste0(OUT, "/Seurat_expr.Robj"))
write.table(x = cellMetaDatax, file = paste0(OUT, "/Seurat_metaData.txt"), row.names = F, col.names = T, quote = F,sep = "\t")
write.table(x = expr@dr$umap@cell.embeddings, file = paste0(OUT, "/Seurat_UMAP_embeddings.txt"), row.names = T, col.names = NA, quote = F,sep = "\t")
#write.table(x = expr.markers_pn, file = paste0(OUT, "/Seurat_markerGenes_pn.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#write.table(x = expr.markers_ftd_labelled, file = paste0(OUT, "/Seurat_markerGenes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#data.table::fwrite(x = as.data.frame(expr@scale.data), file = paste0(OUT, "/UMIcount_scaled.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)

save.image(file = paste0(OUT, "/clustering.RData"))
