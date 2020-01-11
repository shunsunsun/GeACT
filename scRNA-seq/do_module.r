# gene module
setwd("~/lustre/06-Human_cell_atlas/pooled_data/01_stomach/")

library("pheatmap")
library("ggplot2")
library("cowplot")
library("reshape2")
library("dynamicTreeCut")
library("gridExtra")
library("parallel")
library("clusterProfiler")
library("STRINGdb")
source("../../scripts/module_tools.r")
source("../../scripts/pheatmap_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/geneModule/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/module.RData"))

# 1. preprocess ----
# load the dataset
expr_data <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F, comment.char = "")
dim(expr_data)
cellMetaData <- read.table(file = paste0("03-expression/merged/cellCluster/", samplingPos, "/Seurat_metaData.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
dim(cellMetaData)
cellMetaData <- cellMetaData[match(colnames(expr_data), rownames(cellMetaData)), ]

# split the gene expression table by cell type
res <- do_createDT(expr_data, cellMetaData, do.norm = T, cell_num_cutoff = 500)
names(res)

# correlation
res <- do_cor(res, subsp = 500, expr_cutoff = 0.1, mask = F, rm_outlier = F, method = "spearman")

# clustering
res <- do_hc(res, use.abs = F, method = "average")

# region <- grep("^RPS", rownames(res[[1]][["cor_cld"]]))
# pheatmap_new(res[[1]][["cor_cld"]][region, region], cluster_rows = F, cluster_cols = F, 
#              show_rownames = T, show_colnames = F, 
#              breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1), 
#              color = c("blue",colorRampPalette(c("blue", "white", "red"))(100), "red"), gaps_row = NULL, gaps_col = NULL, 
#              highlights_row = c(1,10), highlights_col = c(1,10), highlights_color = "black", display_numbers = F)

# 2. module detection ----
res <- do_detectModule(res, method = "dynamic", h_cutoff = 0.99, size_cutoff = 10)

# do_plotCorHeatmap(res, ctype1 = names(res)[1], ctype2 = names(res)[2], mdid = res[[1]]$cl_list_ftd$cluster[121:125], do_highlights = T, show_rownames = T, show_colnames = F, fontsize_col = 6, vmin = -0.6, vmax = 0.6)
# do_plotHclust(res, ctype1 = 1, mdid = head(res[[1]]$cl_list_ftd$cluster, 5))

# module stat
md_stat <- t(sapply(res, function(x) {
  md_num <- length(unique(x[["cl_table_ftd"]]$cluster))
  gene_num <- nrow(x[["cl_table_ftd"]])
  return(c(md_num, gene_num))
}))
md_stat <- as.data.frame(md_stat)
colnames(md_stat) <- c("module", "gene")
md_stat$ctype <- rownames(md_stat)
md_stat$cellNum <- table(cellMetaData$ident)[rownames(md_stat)]

cl_table_DF <- do.call("rbind", lapply(names(res), function(x) {
  y <- res[[x]][["cl_table_ftd"]][, c("cluster", "size")]
  y <- data.frame(ctype =x , unique(y), stringsAsFactors = F)
}))

pdf(paste0(OUT, "/mdStat.pdf"), width = 4, height = 5)
ggplot(md_stat, aes(x = ctype, y = module)) + geom_bar(fill = "dodgerblue", stat = "identity", show.legend = F) + 
  scale_y_continuous(limits = c(0, max(md_stat$module, md_stat$cellNum/10) * 1.05), expand = c(0, 0), sec.axis = sec_axis(~ . * 10, name = "Cell number")) + 
  geom_line(aes(y = cellNum / 10), group = 1, color = "orange", size = 1.2) + 
  geom_point(data = subset(md_stat, cellNum > 500), aes(y = cellNum / 10)) + 
  geom_hline(yintercept = 500/10, linetype = "dashed") + 
  annotate("text", x = which.min(md_stat$cellNum), y = min(md_stat$cellNum) / 10, label= min(md_stat$cellNum), color = "orange", vjust = -0.5) + 
  annotate("text", x = which.max(md_stat$cellNum), y = max(md_stat$cellNum) / 10, label= max(md_stat$cellNum), color = "orange", vjust = -0.5) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(axis.title.y.right = element_text(color = "orange"), axis.text.y.right = element_text(color = "orange")) + 
  ylab("Module number")

ggplot(cl_table_DF, aes(x = ctype, y = log10(size), fill = ctype)) + geom_violin(show.legend = F) + geom_boxplot(width = 0.2, show.legend = F) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab(expression(log[10] ~ "(Module size)"))

dev.off()

#saveRDS(object = res, file = paste0(OUT, "/res_1.rds"))

# 3. Gene set enrichment ----

# TF
res <- do_enrich(res, db = "TF", expr_cutoff = 10, ncpu = 10)

# miRNA
res <- do_enrich(res, db = "miRNA", expr_cutoff = 0.5, ncpu = 10)

# GO
res <- do_enrichGO(res, ncpu = 10)

# KEGG
res <- do_enrichKEGG(res, ncpu = 10)

# PPI
res <- do_enrichPPI(res, ncpu = 10)

# 4. module map within cell types ----

pdf(paste0(OUT, "/module_map_within_cellType.pdf"), width = 9, height = 6.5)
# do_plotModuleMap(res, ctype = 1, min.avgExpr = 5, min.avgCor = 0.1, showModuleID = T)
for(i in names(res)) {
  do_plotModuleMap(res, ctype = i, showModuleID = F)
  do_plotModuleMap(res, ctype = i, sortBy = "avgExpr", showModuleID = F)
  do_plotModuleMap(res, ctype = i, min.avgExpr = 5, min.avgCor = 0.1)
  do_plotModuleMap(res, ctype = i, sortBy = "avgExpr", min.avgExpr = 5, min.avgCor = 0.1)
}
dev.off()

pdf(paste0(OUT, "/module_map_within_cellType_stat.pdf"), width = 6, height = 5.5)
for(i in 1:length(res)) {
  # expr
  for(db in c("TF", "miRNA", "GO", "KEGG", "PPI")) {
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = avgExpr, y = cvExpr)) + 
      geom_point(aes_string(size = "size", color = paste0("enrich_", db)), alpha = 0.6) + #geom_text(aes(label = cluster), nudge_x = 0.008, size = 3) + 
      ggtitle(names(res)[i])
    print(p)
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = avgExpr)) + geom_density(aes_string(fill = paste0("enrich_", db)), alpha = 0.6) + ggtitle(names(res)[i])
    print(p)
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = cvExpr)) + geom_density(aes_string(fill = paste0("enrich_", db)), alpha = 0.6) + ggtitle(names(res)[i])
    print(p)
  }
  # cor
  for(db in c("TF", "miRNA", "GO", "KEGG", "PPI")) {
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = avgCor, y = cvCor)) + 
      geom_point(aes_string(size = "size", color = paste0("enrich_", db)), alpha = 0.6) + #geom_text(aes(label = cluster), nudge_x = 0.008, size = 3) + 
      ggtitle(names(res)[i])
    print(p)
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = avgCor)) + geom_density(aes_string(fill = paste0("enrich_", db)), alpha = 0.6) + ggtitle(names(res)[i])
    print(p)
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = cvCor)) + geom_density(aes_string(fill = paste0("enrich_", db)), alpha = 0.6) + ggtitle(names(res)[i])
    print(p)
  }
  # expr ~ cor
  for(db in c("TF", "miRNA", "GO", "KEGG", "PPI")) {
    p <- ggplot(res[[i]]$cl_list_ftd, aes(x = avgExpr, y = avgCor)) + 
      geom_point(aes_string(size = "size", color = paste0("enrich_", db)), alpha = 0.6) + #geom_text(aes(label = cluster), nudge_x = 0.008, size = 3) + 
      ggtitle(names(res)[i])
    print(p)
  }
}
dev.off()

# case
pdf(paste0(OUT, "/module_case.pdf"), width = 6, height = 5.5)
for(case_id in res[[1]]$cl_list_ftd$cluster[1:5]) {
  do_plotCorHeatmap(res, ctype1 = names(res)[1], ctype2 = names(res)[2], mdid = case_id, show_rownames = T, fontsize_row = 10, method = "ggplot2", main = paste(names(res)[1], case_id))
  do_plotExprHeatmap(res, ctype1 = names(res)[1], ctype2 = names(res)[2], mdid = case_id)
  do_plotEnrich(res, ctype = names(res)[1], mdid = case_id, db = "GO", main = NULL)
  do_plotEnrich(res, ctype = names(res)[1], mdid = case_id, db = "KEGG", main = NULL)
  do_plotEnrich(res, ctype = names(res)[1], mdid = case_id, db = "MSigDB", category = "C3", main = NULL)
}
rm(case_id)
dev.off()

# each module
pdf(paste0(OUT, "/corHeatmap.pdf"), width = 6, height = 6)
for(ctype in names(res)) { 
  cat(">", ctype, "\n")
  for(mdid in res[[ctype]]$cl_list_ftd$cluster) {
    cat(">>", mdid, "\n")
    annoed <- paste(subset(res[[ctype]]$cl_list_ftd, cluster == mdid, c("enrich_GO", "enrich_MSigDB", "enrich_PPI"), drop = T), collapse = " ")
    do_plotCorHeatmap(res_in = res, ctype1 = ctype, ctype2 = NULL, mdid = mdid, do_highlights = T, 
                      show_rownames = T, fontsize_row = 6, main = paste(ctype, mdid, annoed))
  }
}
dev.off()

pdf(paste0(OUT, "/corHeatmap_flank.pdf"), width = 6, height = 6)
for(ctype in names(res)) { 
  cat(">", ctype, "\n")
  for(mdix in 1:length(res[[ctype]]$cl_list_ftd$cluster)) {
    mdid <- res[[ctype]]$cl_list_ftd$cluster[mdix]
    mdis <- mdid
    if(mdix > 1) { mdis <- c(res[[ctype]]$cl_list_ftd$cluster[mdix-1], mdis) }
    if(mdix < length(res[[ctype]]$cl_list_ftd$cluster)) { mdis <- c(mdis, res[[ctype]]$cl_list_ftd$cluster[mdix+1]) }
    cat(">>", mdid, "\n")
    annoed <- paste(subset(res[[ctype]]$cl_list_ftd, cluster == mdid, c("enrich_GO", "enrich_MSigDB", "enrich_PPI"), drop = T), collapse = " ")
    do_plotCorHeatmap(res_in = res, ctype1 = ctype, ctype2 = NULL, mdid = mdis, do_highlights = T, 
                      show_rownames = T, fontsize_row = 6, main = paste(ctype, mdid, annoed))
  }
}
dev.off()

pdf(paste0(OUT, "/hclust.pdf"), width = 6, height = 6)
for(ctype in names(res)) {
  cat(">", ctype, "\n")
  for(mdid in res[[ctype]]$cl_list_ftd$cluster) {
    cat(">>", mdid, "\n")
    annoed <- paste(subset(res[[ctype]]$cl_list_ftd, cluster == mdid, c("enrich_GO", "enrich_MSigDB", "enrich_PPI"), drop = T), collapse = " ")
    do_plotHclust(res_in = res, ctype1 = ctype, mdid = mdid, horiz = F, do_rev = F, main = paste(ctype, mdid, annoed))
  }
}
dev.off()

pdf(paste0(OUT, "/enrich_GO.pdf"), width = 8, height = 4, onefile = T)
for(cell_type in names(res)) { 
  cat(">", cell_type, "\n")
  for(module_id in res[[cell_type]]$cl_list_ftd$cluster) {
    #cat(">>", module_id, "\n")
    do_plotEnrich(res_in = res, db = "GO", ctype = cell_type, mdid = module_id)
  }
}
dev.off()

# 5. module map across cell types ----
# coordinated plot
#do_plotCoordCor(res, ctype1 = 1, ctype2 = 2, do_plot = T, sub_num = 1:100, show_colnames = T)

# differential module
cairo_pdf(paste0(OUT, "/cmpMd.pdf"), width = 5.5, height = 12, onefile = T)
cmpMd_DF1 <- do_cmpMd(res, ctype1 = names(res)[1], ctype2 = names(res)[2], do_plot = T, fontsize_row = 6, vmin = -0.3, vmax = 0.3)
cmpMd_DF2 <- do_cmpMd(res, ctype1 = names(res)[1], ctype2 = names(res)[3], do_plot = T, fontsize_row = 6, vmin = -0.6, vmax = 0.6)
cmpMd_DF3 <- do_cmpMd(res, ctype1 = names(res)[2], ctype2 = names(res)[3], do_plot = T, fontsize_row = 6, vmin = -0.6, vmax = 0.6)
dev.off()

ggplot(cmpMd_DF1, aes(x = cmpMd_DF1[, 1], y = cmpMd_DF1[, 2], size = -log10(qvalue))) + 
  geom_point(aes(color = factor(holdOn))) + geom_abline(slope = 1, linetype = "dashed") + 
  theme(legend.position = c(0.01, 0.8)) + xlab(names(res)[1]) + ylab(names(res)[2]) + labs(color = "Conserved", size = expression(-log[10] ~ "q-value"))

# plot cor heatmap
# do_plotCorHeatmap(res_in = res, ctype1 = 1, ctype2 = 2, mdid = c("M55","M105","M97"), do_highlights = T, 
#                   fontsize_row = 6, fontsize_col = 6, show_rownames = T)
# do_plotHclust(res, ctype1 = 1, mdid = c("M55","M105","M97"), horiz = F, do_rev = F); abline(h = 1, lty = 2)

# module map
cl_rmdup_LS <- do_mergeModule(res_in = res, ov_cutoff = 0.6, verbose = T)
length(cl_rmdup_LS)
md_map <- t(sapply(cl_rmdup_LS, function(gs) {
  ot <- sapply(names(res), function(x) {
    gsx <- intersect(gs, rownames(res[[x]]$cor_cld))
    tmp <- res[[x]]$cor_cld[gsx, gsx]
    mean_cor <- mean(tmp[upper.tri(tmp)])
    return(mean_cor)
  })
  return(ot)
}))
dim(md_map)
write.table(x = md_map, file = paste0(OUT, "/module_map.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

# attr
md_map_size <- data.frame(mdid= rownames(md_map), size = sapply(rownames(md_map), function(x) { length(cl_rmdup_LS[[x]]) }))
md_map_size$mdid <- factor(md_map_size$mdid, levels = rev(unique(md_map_size$mdid)))
md_map_sorted <- md_map[order(apply(md_map>0.1, 1, sum), decreasing = T), ]

pdf(file = paste0(OUT, "/module_map.pdf"), width = 4, height = 8, onefile = T)
pheatmap_new(md_map_sorted, show_rownames = T, fontsize_row = 2, 
             breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1), 
             color = c("blue",colorRampPalette(c("blue", "white", "red"))(100), "red"), 
             cluster_rows = F, cluster_cols = F , angle_col = 45, cellwidth = 8)
ggplot(md_map_size, aes(x = mdid, y = log10(size), color = "1")) + geom_point(show.legend = F) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 3) + 
  scale_x_discrete(expand = c(0, 2)) + 
  scale_y_continuous(expand = c(0, 0.5)) + coord_flip(clip = "off") + 
  ylab(expression(log[10] ~ "(gene number)")) + 
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 0.5)
dev.off()

pdf(file = paste0(OUT, "/module_map_case.pdf"), width = 5.5, height = 4, onefile = T)
# 1
case_id <- rownames(md_map_sorted)[1]
ggplot(melt(md_map[case_id, , drop = F]), aes(x = Var2, y = Var1, fill = value)) + geom_tile() + geom_text(aes(label = round(value, 2))) + 
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), limits = c(-0.6, 0.6), breaks = c(-0.6, -0.6 / 2, 0, 0.6 / 2, 0.6), na.value = "DDDDDD") + 
  labs(fill = "Correlation") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank()) + 
  theme(legend.margin = margin(l = -5), aspect.ratio = 0.1)
for(nm in seq_along(names(res))) {
  p <- do_plotCorHeatmap(res_in = res, ctype1 = names(res)[nm], ctype2 = NULL, mp_in = cl_rmdup_LS, mpid = case_id, show_rownames = (nm == 1), fontsize_row = 10, 
                         rm.upper = T, method = "ggplot2", main = NULL, show.legend = F) + ylab(names(res)[nm])
  print(p)
}
# 2
case_id <- rownames(md_map_sorted)[2]
ggplot(melt(md_map[case_id, , drop = F]), aes(x = Var2, y = Var1, fill = value)) + geom_tile() + geom_text(aes(label = round(value, 2))) + 
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), limits = c(-0.6, 0.6), breaks = c(-0.6, -0.6 / 2, 0, 0.6 / 2, 0.6), na.value = "DDDDDD") + 
  labs(fill = "Correlation") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank()) + 
  theme(legend.margin = margin(l = -5), aspect.ratio = 0.1)
for(nm in seq_along(names(res))) {
  p <- do_plotCorHeatmap(res_in = res, ctype1 = names(res)[nm], ctype2 = NULL, mp_in = cl_rmdup_LS, mpid = case_id, show_rownames = (nm == 1), fontsize_row = 10, 
                         rm.upper = T, method = "ggplot2", main = NULL, show.legend = F) + ylab(names(res)[nm])
  print(p)
}
# 3
case_id <- rownames(md_map_sorted)[3]
ggplot(melt(md_map[case_id, , drop = F]), aes(x = Var2, y = Var1, fill = value)) + geom_tile() + geom_text(aes(label = round(value, 2))) + 
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), limits = c(-0.6, 0.6), breaks = c(-0.6, -0.6 / 2, 0, 0.6 / 2, 0.6), na.value = "DDDDDD") + 
  labs(fill = "Correlation") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank()) + 
  theme(legend.margin = margin(l = -5), aspect.ratio = 0.1)
for(nm in seq_along(names(res))) {
  p <- do_plotCorHeatmap(res_in = res, ctype1 = names(res)[nm], ctype2 = NULL, mp_in = cl_rmdup_LS, mpid = case_id, show_rownames = (nm == 1), fontsize_row = 10, 
                         rm.upper = T, method = "ggplot2", main = NULL, show.legend = F) + ylab(names(res)[nm])
  print(p)
}
dev.off()

# X. save object ----
save.image(file = paste0(OUT, "/module.RData"))
