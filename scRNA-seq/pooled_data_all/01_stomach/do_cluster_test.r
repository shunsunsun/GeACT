# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/01_stomach/")

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
suppressMessages(library("arrow"))
source("../../scripts/cluster_tools.r")
source("../../scripts/network_tools.r")
source("../../scripts/pheatmap_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/cellCluster/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/clustering.RData"))

load(paste0(OUT, "/Seurat_step2.RData"))

### 3.1 Assigning cell type identity to clusters (ident) ----
cat("> Step 3.1 Running...\n")
expr_assigned <- expr
all(names(expr_assigned@ident) == cellMetaData$cell)
expr_assigned@ident <- factor(paste(cellMetaData$ident, cellMetaData$stage, sep = "."), levels = ct_all)
names(expr_assigned@ident) <- cellMetaData$cell

table(expr_assigned@ident)
#TSNEPlot(object = expr_assigned, do.label = TRUE, pt.size = 1, no.legend = T)

# cell num stat
cellNumStat <- as.data.frame(table(expr_assigned@ident, useNA = "ifany"))
colnames(cellNumStat) <- c("type_stage", "num")
cellNumStat$type <- factor(gsub("\\..*", "", cellNumStat$type_stage), levels = cellTypeMeta_sub$ident)
cellNumStat$stage <- gsub(".*\\.", "", cellNumStat$type_stage)
cellNumStat_cmp <- as.data.frame(acast(cellNumStat, type ~ stage, value.var = "num"))
cellNumStat_cmp$used <- rowSums(cellNumStat_cmp >= 10) == 2
ct_sub <- rownames(cellNumStat_cmp)[cellNumStat_cmp$used]
id_DF <- data.frame(old = ct_sub, new = ct_sub, stringsAsFactors = F)

cl <- makeCluster(5, type = "FORK")
expr.markers_byIdent_LS <- parLapply(cl, ct_sub, function(x) {
  cat(">", x, "\n")
  ct_1 <- paste(x, "19-22w", sep = ".")
  ct_2 <- paste(x, "11-14w", sep = ".")
  expr.markers <- FindMarkers(object = expr_assigned, ident.1 = ct_1, ident.2 = ct_2, test.use = "wilcox", only.pos = F, logfc.threshold = 0, min.pct = 0.25)
  expr.markers$gene <- rownames(expr.markers)
  rownames(expr.markers) <- NULL
  expr.markers$filter <- expr.markers$p_val_adj < 0.05 & abs(expr.markers$avg_logFC) > log(1.5)
  expr.markers$chrXY <- expr.markers$gene %in% c(gene_chrX, gene_chrY)
  expr.markers$cluster <- x
  return(expr.markers)
})
stopCluster(cl); rm(cl)
expr.markers_byIdent_DF <- do.call("rbind", expr.markers_byIdent_LS)
write.table(x = expr.markers_byIdent_DF, file = paste0(OUT, "/Seurat_expr.markers_byIdent.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# plot volcano
pdf(file = file.path(OUT, "DEG_volcano_byIdent.pdf"), width = 4, height = 4)
gp_LS <- lapply(ct_sub, function(x) {
  cat(">", x, "\n")
  expr.markers <- subset(expr.markers_byIdent_DF, cluster == x & ! chrXY)
  expr.markers_sub <- subset(expr.markers, filter)
  expr.markers_sub <- expr.markers_sub[order(expr.markers_sub$avg_logFC), ]
  gp <- ggplot(expr.markers, aes(x = avg_logFC, y = -log10(p_val_adj))) + geom_point(aes(color = filter), alpha = 0.6, show.legend = F) + 
    geom_vline(xintercept = c(-log(1.5), log(1.5)), linetype = "dashed") + 
    coord_cartesian(clip = "off") + 
    scale_x_continuous(limits = c(- max(abs(range(expr.markers$avg_logFC))), max(abs(range(expr.markers$avg_logFC))))) + 
    scale_color_manual(values = c("grey50", "tomato")) + 
    theme(aspect.ratio = 1) + 
    xlab("Log (Fold change)") + ylab(parse(text = "-Log[10]~(FDR)")) + ggtitle(x)
  if(nrow(expr.markers_sub) > 0) {
    idx_1 <- intersect(1:2, 1:nrow(expr.markers_sub))
    idx_2 <- intersect((nrow(expr.markers_sub) - 1):nrow(expr.markers_sub), 1:nrow(expr.markers_sub))
    gp <- gp + geom_text(data = expr.markers_sub[unique(c(idx_1, idx_2)), ], aes(label = gene), nudge_y = -0.5)
  }
  print(gp)
  return(gp)
})
dev.off()

# pdf(file = paste0(OUT, "/Seurat_geneder_diff_gene.pdf"), width = 5, height = 5)
# FeaturePlot(object = expr_assigned, features.plot = c("FOSB"), cols.use = c("grey", "blue"), cells.use = WhichCells(expr_assigned, ident = c(ct_1, ct_2)), reduction.use = "tsne", pt.size = 2, no.legend = T)
# FeaturePlot(object = expr_assigned, features.plot = c("RPS4Y1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# dev.off()

# filter
expr.markers_ftd <- subset(expr.markers_byIdent_DF, filter & ! chrXY)
table(expr.markers_ftd$cluster)
# write
write.table(x = expr.markers_ftd, file = paste0(OUT, "/Seurat_expr.markers_ftd_byIdent.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#
expr.markers_ftd_up <- subset(expr.markers_ftd, avg_logFC > 0)
expr.markers_ftd_dw <- subset(expr.markers_ftd, avg_logFC < 0)
#
top10_up <- expr.markers_ftd_up %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10_dw <- expr.markers_ftd_dw %>% group_by(cluster) %>% top_n(10, -avg_logFC)

# GO enrichment for obtained gene markers
expr.markers_ftd_updw <- expr.markers_ftd
expr.markers_ftd_updw$cluster <- paste(expr.markers_ftd_updw$cluster, ifelse(expr.markers_ftd_updw$avg_logFC > 0, "up", "down"), sep = ".")

enriched_byIdent_LS <- do_GOenrich(expr.markers_ftd_updw, ncpu = 8) # slow

topEnrich <- do.call("rbind", lapply(enriched_byIdent_LS, function(x) { y <- x[1:20, c(9,1:8)] }))
#View(topEnrich)

# plot GO enrich
enrich_res <- do.call("rbind", enriched_byIdent_LS)
enrich_res$Term <- Hmisc::capitalize(gsub(" \\(.*\\)", "", enrich_res$Term))
enrich_res$Adjusted.P.value <- enrich_res$q_value
#enrich_res <- enrich_res[order(enrich_res$Adjusted.P.value), ]
enrich_res <- subset(enrich_res, Adjusted.P.value < 0.05)
enrich_res$type <- gsub(".*\\.", "", enrich_res$cluster)
enrich_res$cluster <- gsub("\\..*", "", enrich_res$cluster)
# write
write.table(x = enrich_res, file = file.path(OUT, "DEG_GOenrich_byIdent.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#

# 1) up
enrich_res_sub <- subset(enrich_res, type == "up")
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res_sub, enrich_res_sub$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, id_DF, by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# add dummy
nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new),
                                                 Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
enrich_res_DF$new <- factor(enrich_res_DF$new, levels = id_DF$new)
enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", samplingPos, "/DEG_up_GOenrich_byIdent.pdf"), width = 12, height = 8)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) +
  geom_bar(stat = "identity", width = 0.9, show.legend = F) +
  facet_grid(new ~ ., scales ="free", drop = F) +
  coord_flip() +
  ylab(expression(paste(-Log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + ggtitle("Up regulated genes") + 
  theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.01)) +
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") +
  theme(panel.spacing.y = unit(0.2, "lines"))
gp
dev.off()

# 2) down
enrich_res_sub <- subset(enrich_res, type == "down")
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res_sub, enrich_res_sub$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, id_DF, by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# add dummy
nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new),
                                                 Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
enrich_res_DF$new <- factor(enrich_res_DF$new, levels = id_DF$new)
enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", samplingPos, "/DEG_down_GOenrich_byIdent.pdf"), width = 12, height = 8)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) +
  geom_bar(stat = "identity", width = 0.9, show.legend = F) +
  facet_grid(new ~ ., scales ="free", drop = F) +
  coord_flip() +
  ylab(expression(paste(-Log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + ggtitle("Down regulated genes") + 
  theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.01)) +
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") +
  theme(panel.spacing.y = unit(0.2, "lines"))
gp
dev.off()

# point
# tmp <- expr_data_avg
# tmp$DEG <- rownames(tmp) %in% expr.markers_ftd$gene
# gp <- ggplot(tmp, aes(x = log10(`14w` + 1), y = log10(`20w` + 1), color = DEG)) + geom_point(alpha = 0.6, show.legend = F) + 
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
#   scale_color_manual(values = c("grey60", "dodgerblue3")) + 
#   xlab("~14-week-old") + ylab("~20-week-old") + 
#   theme(aspect.ratio = 1)
# ggsave(filename = paste0(OUT, "/DEG_point_14w_vs_20w.pdf"), plot = gp, width = 5, height = 5)

# # heatmap
# expr.markers_ftd_labelled <- expr.markers_ftd
# expr.markers_ftd_labelled$ident <- factor(expr.markers_ftd_labelled$cluster, levels = ct_sub)
# tmp <- expr.markers_ftd_labelled %>% group_by(cluster) %>% top_n(10, avg_logFC)
# tmp <- tmp[! duplicated(tmp$gene), ]
# expr_assigned_tmp <- expr_assigned
# #expr_assigned_tmp@ident <- factor(cellMetaData$stage)
# #names(expr_assigned_tmp@ident) <- names(expr_assigned@ident)
# gp <- DoHeatmap_new(object = expr_assigned_tmp, genes.use = tmp$gene, genes.group = tmp$ident, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
#               group.cex = 11, cex.row = 7, do.colBar = F, strip.text.x.top = 15, group.label.rot = T, 
#               panel.spacing.y = 0, 
#               cells.use = WhichCells(expr_assigned, max.cells.per.ident = 200, random.seed = 1))
# ggsave(filename = paste0(OUT, "/DEG_heatmap_byIdent.pdf"), plot = gp, width = 16, height = 8)

# TF
library("RcisTarget")

geneLists <- split(expr.markers_ftd_updw$gene, expr.markers_ftd_updw$cluster)
motifRankings <- importRankings("~/lustre/06-Human_cell_atlas/Data/RcisTarget/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
data(motifAnnotations_hgnc)
###
motifEnrichmentTable_wGenes_byIdent <- cisTarget(geneLists, motifRankings, motifAnnot = motifAnnotations_hgnc) # slow
###
motifEnrichmentTable_wGenes <- as.data.frame(motifEnrichmentTable_wGenes_byIdent)
# write
write.table(x = motifEnrichmentTable_wGenes, file = file.path(OUT, "DEG_motifEnrich_byIdent.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#
#motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
#DT::datatable(motifEnrichmentTable_wGenes_wLogo[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], escape = FALSE, filter="top", options = list(pageLength = 5))

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf, motifEnrichmentTable_wGenes$geneSet), function(x) {
  genes <- gsub(" \\(.*\\).", "; ", x, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  return(genesSplit)
})
#anotatedTfs

network_LS <- apply(motifEnrichmentTable_wGenes, 1, function(x) {
  TF <- gsub(" \\(.*\\).", "", x[5], fixed=FALSE)
  TF <- unlist(strsplit(TF, "; "))
  target <- unlist(strsplit(x[9], split = ";"))
  network <- data.table::CJ(TF, target)
  network <- data.frame(geneSet = rep(x[1], nrow(network)), motif = rep(x[2], nrow(network)), 
                        NES = as.numeric(rep(x[3], nrow(network))), AUC = as.numeric(rep(x[4], nrow(network))), 
                        network, stringsAsFactors = F)
  return(network)
})
network_DF <- do.call("rbind", network_LS)

# need TF differentially expressed
network_ftd_LS <- lapply(unique(network_DF$geneSet), function(x) {
  network_DF_sub <- subset(network_DF, geneSet == x)
  x_cluster <- gsub("\\..*", "", x)
  expr.markers_ftd_sub <- subset(expr.markers_ftd, cluster == x_cluster)
  y <- subset(network_DF_sub, TF %in% expr.markers_ftd_sub$gene)
  #y$TF.avg_logFC <- expr.markers_ftd_sub$avg_logFC[match(y$TF, expr.markers_ftd_sub$gene)]
  #y$TF.p_val_adj <- expr.markers_ftd_sub$p_val_adj[match(y$TF, expr.markers_ftd_sub$gene)]
  return(y)
})
network_ftd_DF <- do.call("rbind", network_ftd_LS)
table(network_ftd_DF$geneSet)

# plot TF expression
# VlnPlot(object = expr_assigned_tmp, features.plot = expr.markers_ftd$gene[1:4], nCol = 2, point.size.use = 0.01)

# plot network
suppressMessages(library("igraph"))

pdf(file = file.path(OUT, "DEG_TF_network_byIdent.pdf"), width = 6, height = 6)

for(i in unique(gsub("\\..*", "", network_ftd_DF$geneSet))) {
  do_plotNt(gs_case = i, show.all.nodes = F, print.nodes = F)
}

dev.off()
##

### 3.2 Assigning cell type identity to clusters (group) ----
cat("> Step 3.2 Running...\n")
expr_assigned <- expr
all(names(expr_assigned@ident) == cellMetaData$cell)
expr_assigned@ident <- factor(paste(cellMetaData$group, cellMetaData$stage, sep = "."), levels = cg_all)
names(expr_assigned@ident) <- cellMetaData$cell

table(expr_assigned@ident)
#TSNEPlot(object = expr_assigned, do.label = TRUE, pt.size = 1, no.legend = T)

# cell num stat
cellNumStat <- as.data.frame(table(expr_assigned@ident, useNA = "ifany"))
colnames(cellNumStat) <- c("type_stage", "num")
cellNumStat$type <- factor(gsub("\\..*", "", cellNumStat$type_stage), levels = unique(cellTypeMeta_sub$group))
cellNumStat$stage <- gsub(".*\\.", "", cellNumStat$type_stage)
cellNumStat_cmp <- as.data.frame(acast(cellNumStat, type ~ stage, value.var = "num"))
cellNumStat_cmp$used <- rowSums(cellNumStat_cmp >= 10) == 2
ct_sub <- rownames(cellNumStat_cmp)[cellNumStat_cmp$used]
id_DF <- data.frame(old = ct_sub, new = ct_sub, stringsAsFactors = F)

# DEGs
cl <- makeCluster(5, type = "FORK")
expr.markers_byGroup_LS <- parLapply(cl, ct_sub, function(x) {
  cat(">", x, "\n")
  ct_1 <- paste(x, "19-22w", sep = ".")
  ct_2 <- paste(x, "11-14w", sep = ".")
  expr.markers <- FindMarkers(object = expr_assigned, ident.1 = ct_1, ident.2 = ct_2, test.use = "wilcox", only.pos = F, logfc.threshold = 0, min.pct = 0.25)
  expr.markers$gene <- rownames(expr.markers)
  rownames(expr.markers) <- NULL
  expr.markers$filter <- expr.markers$p_val_adj < 0.05 & abs(expr.markers$avg_logFC) > log(1.5)
  expr.markers$chrXY <- expr.markers$gene %in% c(gene_chrX, gene_chrY)
  expr.markers$cluster <- x
  return(expr.markers)
})
stopCluster(cl); rm(cl)
expr.markers_byGroup_DF <- do.call("rbind", expr.markers_byGroup_LS)
write.table(x = expr.markers_byGroup_DF, file = paste0(OUT, "/Seurat_expr.markers_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# plot volcano
pdf(file = file.path(OUT, "DEG_volcano_byGroup.pdf"), width = 4, height = 4)
gp_LS <- lapply(ct_sub, function(x) {
  cat(">", x, "\n")
  expr.markers <- subset(expr.markers_byGroup_DF, cluster == x & ! chrXY)
  expr.markers_sub <- subset(expr.markers, filter)
  expr.markers_sub <- expr.markers_sub[order(expr.markers_sub$avg_logFC), ]
  gp <- ggplot(expr.markers, aes(x = avg_logFC, y = -log10(p_val_adj))) + geom_point(aes(color = filter), alpha = 0.6, show.legend = F) + 
    geom_vline(xintercept = c(-log(1.5), log(1.5)), linetype = "dashed") + 
    coord_cartesian(clip = "off") + 
    scale_x_continuous(limits = c(- max(abs(range(expr.markers$avg_logFC))), max(abs(range(expr.markers$avg_logFC))))) + 
    scale_color_manual(values = c("grey50", "tomato")) + 
    theme(aspect.ratio = 1) + 
    xlab("Log (Fold change)") + ylab(parse(text = "-Log[10]~(FDR)")) + ggtitle(x)
  if(nrow(expr.markers_sub) > 0) {
    idx_1 <- intersect(1:2, 1:nrow(expr.markers_sub))
    idx_2 <- intersect((nrow(expr.markers_sub) - 1):nrow(expr.markers_sub), 1:nrow(expr.markers_sub))
    gp <- gp + geom_text(data = expr.markers_sub[unique(c(idx_1, idx_2)), ], aes(label = gene), nudge_y = -0.5)
  }
  print(gp)
  return(gp)
})
dev.off()

# pdf(file = paste0(OUT, "/Seurat_geneder_diff_gene.pdf"), width = 5, height = 5)
# FeaturePlot(object = expr_assigned, features.plot = c("FOSB"), cols.use = c("grey", "blue"), cells.use = WhichCells(expr_assigned, ident = c(ct_1, ct_2)), reduction.use = "tsne", pt.size = 2, no.legend = T)
# FeaturePlot(object = expr_assigned, features.plot = c("RPS4Y1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# dev.off()

# filter
expr.markers_ftd <- subset(expr.markers_byGroup_DF, filter & ! chrXY)
table(expr.markers_ftd$cluster)
# write
write.table(x = expr.markers_ftd, file = paste0(OUT, "/Seurat_expr.markers_ftd_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#
expr.markers_ftd_up <- subset(expr.markers_ftd, avg_logFC > 0)
expr.markers_ftd_dw <- subset(expr.markers_ftd, avg_logFC < 0)
#
top10_up <- expr.markers_ftd_up %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10_dw <- expr.markers_ftd_dw %>% group_by(cluster) %>% top_n(10, -avg_logFC)

# GO enrichment for obtained gene markers
expr.markers_ftd_updw <- expr.markers_ftd
expr.markers_ftd_updw$cluster <- paste(expr.markers_ftd_updw$cluster, ifelse(expr.markers_ftd_updw$avg_logFC > 0, "up", "down"), sep = ".")

enriched_byGroup_LS <- do_GOenrich(expr.markers_ftd_updw, ncpu = 8) # slow

topEnrich <- do.call("rbind", lapply(enriched_byGroup_LS, function(x) { y <- x[1:20, c(9,1:8)] }))
#View(topEnrich)

# plot GO enrich
enrich_res <- do.call("rbind", enriched_byGroup_LS)
enrich_res$Term <- Hmisc::capitalize(gsub(" \\(.*\\)", "", enrich_res$Term))
enrich_res$Adjusted.P.value <- enrich_res$q_value
#enrich_res <- enrich_res[order(enrich_res$Adjusted.P.value), ]
enrich_res <- subset(enrich_res, Adjusted.P.value < 0.05)
enrich_res$type <- gsub(".*\\.", "", enrich_res$cluster)
enrich_res$cluster <- gsub("\\..*", "", enrich_res$cluster)
# write
write.table(x = enrich_res, file = file.path(OUT, "DEG_GOenrich_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#

# 1) up
enrich_res_sub <- subset(enrich_res, type == "up")
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res_sub, enrich_res_sub$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, id_DF, by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# add dummy
nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new),
                                                 Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
enrich_res_DF$new <- factor(enrich_res_DF$new, levels = id_DF$new)
enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", samplingPos, "/DEG_up_GOenrich_byGroup.pdf"), width = 12, height = 8)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) +
  geom_bar(stat = "identity", width = 0.9, show.legend = F) +
  facet_grid(new ~ ., scales ="free", drop = F) +
  coord_flip() +
  ylab(expression(paste(-Log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + ggtitle("Up regulated genes") + 
  theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.01)) +
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") +
  theme(panel.spacing.y = unit(0.2, "lines"))
gp
dev.off()

# 2) down
enrich_res_sub <- subset(enrich_res, type == "down")
enrich_res_DF <- do.call("rbind", lapply(split(enrich_res_sub, enrich_res_sub$cluster), function(x) { head(x, 2) }))
enrich_res_DF <- merge(enrich_res_DF, id_DF, by.x = "cluster", by.y = "old", sort = F)
enrich_res_DF <- enrich_res_DF[, c("new", "Term", "Adjusted.P.value")]

# add dummy
nsnum <- length(setdiff(id_DF$new, enrich_res_DF$new))
enrich_res_DF <- rbind(enrich_res_DF, data.frame(new=setdiff(id_DF$new, enrich_res_DF$new),
                                                 Term=rep("N.S.", nsnum), Adjusted.P.value = rep(1, nsnum), stringsAsFactors = F))
# sort cell type
enrich_res_DF$new <- factor(enrich_res_DF$new, levels = id_DF$new)
enrich_res_DF$Term <- factor(enrich_res_DF$Term, levels = rev(unique(enrich_res_DF$Term)))

pdf(file = paste0(OUT, "/", samplingPos, "/DEG_down_GOenrich_byGroup.pdf"), width = 12, height = 8)
gp <- ggplot(enrich_res_DF, aes(x = Term, y = -log10(Adjusted.P.value), fill = new)) +
  geom_bar(stat = "identity", width = 0.9, show.legend = F) +
  facet_grid(new ~ ., scales ="free", drop = F) +
  coord_flip() +
  ylab(expression(paste(-Log[10], " (Adjusted P)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + ggtitle("Down regulated genes") + 
  theme(strip.text.y = element_text(angle = 0, margin = margin(l = 5, r = 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(enrich_res_DF$Adjusted.P.value)) * 1.01)) +
  theme(axis.ticks.y = element_blank(), strip.background = element_blank(), strip.text.y = element_text(hjust = 0)) + xlab("GO term") +
  theme(panel.spacing.y = unit(0.2, "lines"))
gp
dev.off()

# point
# tmp <- expr_data_avg
# tmp$DEG <- rownames(tmp) %in% expr.markers_ftd$gene
# gp <- ggplot(tmp, aes(x = log10(`14w` + 1), y = log10(`20w` + 1), color = DEG)) + geom_point(alpha = 0.6, show.legend = F) + 
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
#   scale_color_manual(values = c("grey60", "dodgerblue3")) + 
#   xlab("~14-week-old") + ylab("~20-week-old") + 
#   theme(aspect.ratio = 1)
# ggsave(filename = paste0(OUT, "/DEG_point_14w_vs_20w.pdf"), plot = gp, width = 5, height = 5)

# # heatmap
# expr.markers_ftd_labelled <- expr.markers_ftd
# expr.markers_ftd_labelled$ident <- factor(expr.markers_ftd_labelled$cluster, levels = ct_sub)
# tmp <- expr.markers_ftd_labelled %>% group_by(cluster) %>% top_n(10, avg_logFC)
# tmp <- tmp[! duplicated(tmp$gene), ]
# expr_assigned_tmp <- expr_assigned
# #expr_assigned_tmp@ident <- factor(cellMetaData$stage)
# #names(expr_assigned_tmp@ident) <- names(expr_assigned@ident)
# gp <- DoHeatmap_new(object = expr_assigned_tmp, genes.use = tmp$gene, genes.group = tmp$ident, slim.col.label = TRUE, remove.key = F, rotate.key = T, 
#               group.cex = 11, cex.row = 7, do.colBar = F, strip.text.x.top = 15, group.label.rot = T, 
#               panel.spacing.y = 0, 
#               cells.use = WhichCells(expr_assigned, max.cells.per.ident = 200, random.seed = 1))
# ggsave(filename = paste0(OUT, "/DEG_heatmap_byGroup.pdf"), plot = gp, width = 16, height = 8)

# TF
library("RcisTarget")

geneLists <- split(expr.markers_ftd_updw$gene, expr.markers_ftd_updw$cluster)
motifRankings <- importRankings("~/lustre/06-Human_cell_atlas/Data/RcisTarget/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
data(motifAnnotations_hgnc)
###
motifEnrichmentTable_wGenes_byGroup <- cisTarget(geneLists, motifRankings, motifAnnot = motifAnnotations_hgnc) # slow
###
motifEnrichmentTable_wGenes <- as.data.frame(motifEnrichmentTable_wGenes_byGroup)
# write
write.table(x = motifEnrichmentTable_wGenes, file = file.path(OUT, "DEG_motifEnrich_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#
#motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
#DT::datatable(motifEnrichmentTable_wGenes_wLogo[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], escape = FALSE, filter="top", options = list(pageLength = 5))

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf, motifEnrichmentTable_wGenes$geneSet), function(x) {
  genes <- gsub(" \\(.*\\).", "; ", x, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  return(genesSplit)
})
#anotatedTfs

network_LS <- apply(motifEnrichmentTable_wGenes, 1, function(x) {
  TF <- gsub(" \\(.*\\).", "", x[5], fixed=FALSE)
  TF <- unlist(strsplit(TF, "; "))
  target <- unlist(strsplit(x[9], split = ";"))
  network <- data.table::CJ(TF, target)
  network <- data.frame(geneSet = rep(x[1], nrow(network)), motif = rep(x[2], nrow(network)), 
                        NES = as.numeric(rep(x[3], nrow(network))), AUC = as.numeric(rep(x[4], nrow(network))), 
                        network, stringsAsFactors = F)
  return(network)
})
network_DF <- do.call("rbind", network_LS)

# need TF differentially expressed
network_ftd_LS <- lapply(unique(network_DF$geneSet), function(x) {
  network_DF_sub <- subset(network_DF, geneSet == x)
  x_cluster <- gsub("\\..*", "", x)
  expr.markers_ftd_sub <- subset(expr.markers_ftd, cluster == x_cluster)
  y <- subset(network_DF_sub, TF %in% expr.markers_ftd_sub$gene)
  #y$TF.avg_logFC <- expr.markers_ftd_sub$avg_logFC[match(y$TF, expr.markers_ftd_sub$gene)]
  #y$TF.p_val_adj <- expr.markers_ftd_sub$p_val_adj[match(y$TF, expr.markers_ftd_sub$gene)]
  return(y)
})
network_ftd_DF <- do.call("rbind", network_ftd_LS)
table(network_ftd_DF$geneSet)

# plot TF expression
# VlnPlot(object = expr_assigned_tmp, features.plot = expr.markers_ftd$gene[1:4], nCol = 2, point.size.use = 0.01)

# plot network
suppressMessages(library("igraph"))

pdf(file = file.path(OUT, "DEG_TF_network_byGroup.pdf"), width = 6, height = 6)

for(i in unique(gsub("\\..*", "", network_ftd_DF$geneSet))) {
  do_plotNt(gs_case = i, show.all.nodes = F, print.nodes = F)
}

dev.off()
##

# 4. save and write meta table ----
# save(expr, file = paste0(OUT, "/Seurat_expr.Robj"))
# combine meta with identity
# cellMetaDatax <- cellMetaData
# write.table(x = cellMetaDatax, file = paste0(OUT, "/Seurat_metaData.txt"), row.names = F, col.names = T, quote = F,sep = "\t")
# write.table(x = expr.markers_pn, file = paste0(OUT, "/Seurat_markerGenes_pn.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
# write.table(x = expr.markers_ftd_labelled, file = paste0(OUT, "/Seurat_markerGenes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
# data.table::fwrite(x = as.data.frame(expr@scale.data), file = paste0(OUT, "/UMIcount_scaled.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)

save.image(file = paste0(OUT, "/clustering.RData"))
