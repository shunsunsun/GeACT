# gene module
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

suppressMessages({
  library("pheatmap")
  library("ggplot2")
  library("cowplot")
  library("reshape2")
  library("dynamicTreeCut")
  library("gridExtra")
  library("ggalluvial")
  library("parallel")
  library("clusterProfiler")
  library("STRINGdb")
  library("ComplexHeatmap")
})
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
cellMetaData <- read.table(file = paste0("03-expression/merged/cellCluster/", samplingPos, "/Seurat_metaData_pooled.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
dim(cellMetaData)
cell_metatable <- read.table("cell_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
dim(cell_metatable)
cellMetaData <- merge(cellMetaData, cell_metatable[, "tissue", drop = F], by = 0, sort = F)
cellMetaData$ident.ori <- cellMetaData$ident
cellMetaData$ident <- Hmisc::capitalize(paste(cellMetaData$tissue, cellMetaData$ident, sep = "."))
rownames(cellMetaData) <- cellMetaData$Row.names
cellMetaData <- cellMetaData[, -1]

# split the gene expression table by cell type
res <- do_createDT(expr_data, cellMetaData, do.norm = T, cell_num_cutoff = 500)

# correlation
res <- do_cor(res, subsp = 500, expr_cutoff = 0.1, mask = F, rm_outlier = F, method = "spearman")

# HC clustering
res <- do_hc(res, use.abs = F, method = "average")

# region <- grep("^RPS", rownames(res[[1]][["cor_cld"]]))
# pheatmap_new(res[[1]][["cor_cld"]][region, region], cluster_rows = F, cluster_cols = F, 
#              show_rownames = T, show_colnames = F, 
#              breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1), 
#              color = c("blue",colorRampPalette(c("blue", "white", "red"))(100), "red"), gaps_row = NULL, gaps_col = NULL, 
#              highlights_row = c(1,10), highlights_col = c(1,10), highlights_color = "black", display_numbers = F)

# 2. module detection ----
res <- do_detectModule(res, method = "dynamic", h_cutoff = 0.99, size_cutoff = 10, avgCor_cutoff = 0.088)

# do_plotCorHeatmap(res, ctype1 = 1, ctype2 = 2, mdid = head(res[[1]]$cl_list_ftd$cluster, 5), do_highlights = T, show_colnames = T, fontsize_col = 6, vmin = -0.6, vmax = 0.6)
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
  y <- res[[x]][["cl_table_ftd"]][, c("cluster", "size")]; y <- data.frame(ctype =x , unique(y), stringsAsFactors = F)
}))

pdf(paste0(OUT, "/mdStat.pdf"), width = 6, height = 6)
ggplot(md_stat, aes(x = ctype, y = module)) + geom_bar(fill = "dodgerblue", stat = "identity", show.legend = F) + 
  scale_y_continuous(limits = c(0, max(md_stat$module, md_stat$cellNum/10) * 1.05), expand = c(0, 0), sec.axis = sec_axis(~ . * 10, name = "Cell number")) + 
  geom_line(aes(y = cellNum / 10), group = 1, color = "orange", size = 1.2) + 
  geom_point(data = subset(md_stat, cellNum > 500), aes(y = cellNum / 10)) + 
  #geom_hline(yintercept = 500/10, linetype = "dashed") + 
  annotate("text", x = which.min(md_stat$cellNum), y = min(md_stat$cellNum) / 10, label= min(md_stat$cellNum), color = "orange", vjust = 1.5) + 
  annotate("text", x = which.max(md_stat$cellNum), y = max(md_stat$cellNum) / 10, label= max(md_stat$cellNum), color = "orange", vjust = -0.5) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(axis.title.y.right = element_text(color = "orange"), axis.text.y.right = element_text(color = "orange")) + 
  ylab("Module number")

ggplot(cl_table_DF, aes(x = ctype, y = log10(size), fill = ctype)) + geom_violin(show.legend = F) + geom_boxplot(width = 0.2, show.legend = F) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab(expression(log[10] ~ " Module size"))

dev.off()

#saveRDS(object = res, file = paste0(OUT, "/res_1.rds"))

# 3. Gene set enrichment ----

# TF
# res <- do_enrich(res, db = "TF", expr_cutoff = 10, ncpu = 5)

# miRNA
# res <- do_enrich(res, db = "miRNA", expr_cutoff = 0.5, ncpu = 5)

# GO
# res <- do_enrichGO(res, ncpu = 5)

# KEGG
# res <- do_enrichKEGG(res, ncpu = 5)

# PPI
# res <- do_enrichPPI(res, ncpu = 5)

# 4. module map within cell types ----

pdf(paste0(OUT, "/module_map_within_cellType.pdf"), width = 9, height = 6.5)
# do_plotModuleMap(res, ctype = 1, min.avgExpr = 5, min.avgCor = 0.1, showModuleID = T)
for(i in names(res)) {
  do_plotModuleMap(res, ctype = i, showModuleID = F)
  #do_plotModuleMap(res, ctype = i, sortBy = "avgExpr", showModuleID = F)
  #do_plotModuleMap(res, ctype = i, min.avgExpr = 5, min.avgCor = 0.1)
  #do_plotModuleMap(res, ctype = i, sortBy = "avgExpr", min.avgExpr = 5, min.avgCor = 0.1)
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
# 1
pdf(paste0(OUT, "/module_case_ct1.pdf"), width = 8, height = 5.5)
do_plotCorHeatmap(res, ctype1 = names(res)[1], mdid = "M53", show_rownames = T, fontsize_row = 6)
do_plotExprHeatmap(res, ctype1 = names(res)[1], mdid = "M53", show_rownames = T)
do_plotHclust(res, ctype1 = names(res)[1], mdid = "M53")
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M202", db = "GO", main = NULL)
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M202", db = "KEGG", main = NULL)
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M202", db = "MSigDB", category = "C3", main = NULL)
dev.off()
# 2
pdf(paste0(OUT, "/module_case_ct2.pdf"), width = 8, height = 5.5)
do_plotCorHeatmap(res, ctype1 = names(res)[2], mdid = "M13", show_rownames = T, fontsize_row = 6)
do_plotExprHeatmap(res, ctype1 = names(res)[2], mdid = "M13", show_rownames = T)
do_plotHclust(res, ctype1 = names(res)[2], mdid = "M13")
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M23", db = "GO", charMaxLen = 32, main = NULL)
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M23", db = "KEGG", main = NULL)
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M23", db = "MSigDB", category = "C3", main = NULL)
dev.off()
# 3
pdf(paste0(OUT, "/module_case_ct3.pdf"), width = 8, height = 5.5)
do_plotCorHeatmap(res, ctype1 = names(res)[3], mdid = "M36", show_rownames = T, fontsize_row = 6)
do_plotExprHeatmap(res, ctype1 = names(res)[3], mdid = "M36", show_rownames = T)
do_plotHclust(res, ctype1 = names(res)[3], mdid = "M36")
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M137", db = "GO", charMaxLen = 32, main = NULL)
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M137", db = "KEGG", main = NULL)
#do_plotEnrich(res, ctype = names(res)[1], mdid = "M137", db = "MSigDB", category = "C3", main = NULL)
dev.off()

# each module
pdf(paste0(OUT, "/corHeatmap.pdf"), width = 6, height = 6)
for(ctype in names(res)) { 
  cat(">", ctype, "\n")
  for(mdid in res[[ctype]]$cl_list_ftd$cluster) {
    cat(">>", mdid, "\n")
    annoed <- "X" #paste(subset(res[[ctype]]$cl_list_ftd, cluster == mdid, c("enrich_GO", "enrich_MSigDB", "enrich_PPI"), drop = T), collapse = " ")
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

pdf(paste0(OUT, "/corHeatmap_hclust.pdf"), width = 6, height = 6)
for(ctype in names(res)[1]) { 
  cat(">", ctype, "\n")
  for(mdid in res[[ctype]]$cl_list_ftd$cluster) {
    cat(">>", mdid, "\n")
    annoed <- "X" #paste(subset(res[[ctype]]$cl_list_ftd, cluster == mdid, c("enrich_GO", "enrich_MSigDB", "enrich_PPI"), drop = T), collapse = " ")
    do_plotCorHeatmap(res_in = res, ctype1 = ctype, ctype2 = NULL, mdid = mdid, do_highlights = T, 
                      show_rownames = T, fontsize_row = 6, main = paste(ctype, mdid, annoed))
    do_plotHclust(res_in = res, ctype1 = ctype, mdid = mdid, horiz = F, do_rev = F, main = paste(ctype, mdid, annoed))
    abline(h = 1 - 0.0875, lty = 2)
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
# cairo_pdf(paste0(OUT, "/cmpMd.pdf"), width = 5.5, height = 12, onefile = T)
# cmpMd_DF1 <- do_cmpMd(res, ctype1 = names(res)[1], ctype2 = names(res)[2], do_plot = T, fontsize_row = 6, vmin = -0.3, vmax = 0.3)
# cmpMd_DF2 <- do_cmpMd(res, ctype1 = names(res)[1], ctype2 = names(res)[3], do_plot = T, fontsize_row = 6, vmin = -0.6, vmax = 0.6)
# cmpMd_DF3 <- do_cmpMd(res, ctype1 = names(res)[2], ctype2 = names(res)[3], do_plot = T, fontsize_row = 6, vmin = -0.6, vmax = 0.6)
# dev.off()

# ggplot(cmpMd_DF1, aes(x = cmpMd_DF1$`FB-ADAM28`, y = cmpMd_DF1$`FB-FBLN1`, size = -log10(qvalue))) + 
#   geom_point(aes(color = factor(holdOn))) + geom_abline(slope = 1, linetype = "dashed") + 
#   theme(legend.position = c(0.01, 0.8)) + xlab(names(res)[1]) + ylab(names(res)[2]) + labs(color = "Conserved", size = expression(-log[10] ~ "q-value"))

# plot cor heatmap
do_plotCorHeatmap(res_in = res, ctype1 = 1, ctype2 = 2, mdid = c("M53","M22","M109"), do_highlights = T, 
                  fontsize_row = 6, fontsize_col = 6, show_rownames = T)
do_plotHclust(res, ctype1 = 1, mdid = c("M53","M22","M109"), horiz = F, do_rev = F); abline(h = 1, lty = 2)

# module map
mes <- do_mergeModule(res_in = res, ov_cutoff = 0.85, verbose = T)
cl_rmdup_LS <- split(mes$merged$cl_table_ftd$gene, mes$merged$cl_table_ftd$cluster)
cl_rmdup_LS <- cl_rmdup_LS[unique(mes$merged$cl_table_ftd$cluster)]
length(cl_rmdup_LS)

# write output
write.table(x = mes$merged$cl_table_ftd, file = paste0(OUT, "/geneModule_", "merged", ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
md_avgExpr <- mes$merged$avgExpr
write.table(x = md_avgExpr, file = paste0(OUT, "/module_map_avgExpr.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
md_map <- mes$merged$avgCor
write.table(x = md_map, file = paste0(OUT, "/module_map_avgCor.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

# gene set enrichment
# TF
mes <- do_enrich(mes, db = "TF", expr_cutoff = 10, res_in_sg = res, ncpu = 8)

# miRNA
mes <- do_enrich(mes, db = "miRNA", expr_cutoff = 0.5, res_in_sg = res, ncpu = 8)

# GO
mes <- do_enrichGO(mes, ncpu = 8)

# KEGG
mes <- do_enrichKEGG(mes, ncpu = 8)

# PPI
mes <- do_enrichPPI(mes, ncpu = 8)

mes <- do_enrichEdge(res_in = mes, db = "HIPPIE", category = "PPI", p_cutoff = 1e-10, ncpu = 5)

mes <- do_enrichEdge(res_in = mes, db = "HuRI", category = "PPI", p_cutoff = 0.05, ncpu = 5)

saveRDS(mes, file = paste0(OUT, "/mes.rds"))

# attr
md_map_size <- data.frame(mdid= rownames(md_map), size = sapply(rownames(md_map), function(x) { length(cl_rmdup_LS[[x]]) }))
md_map_size$mdid <- factor(md_map_size$mdid, levels = rev(unique(md_map_size$mdid)))
md_map_sorted <- md_map[order(apply(md_map >= 0.088, 1, sum), decreasing = T), ]

GO_top5 <- lapply(split(mes$merged$enrich_GO, mes$merged$enrich_GO$cluster), function(x) { y <- head(x, 5) })
GO_top5 <- GO_top5[unique(mes$merged$enrich_GO$cluster)]
GO_top5 <- do.call("rbind", GO_top5)
#View(GO_top5)

# simple Heatmap
# pheatmap_new(md_map_sorted, show_rownames = F, fontsize_row = 2, fontsize_col = 8, 
#              breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1), 
#              color = c("blue",colorRampPalette(c("blue", "white", "red"))(100), "red"), 
#              cluster_rows = F, cluster_cols = T , angle_col = 45, cellwidth = 10)

# complex Heatmap
mdmap_ts <- gsub("_", " ", sapply(strsplit(colnames(md_map), split = ".", fixed = T), '[', 1))
mdmap_ct <- gsub("_", " ", sapply(strsplit(colnames(md_map), split = ".", fixed = T), '[', 2))
mdmap_ds <- sapply(split(mes$merged$enrich_GO$Description, mes$merged$enrich_GO$cluster), function(xa) { 
  x <- xa[1:10]
  if(any(grepl("mitochondrial", x))) { y <- "mitochondria" }
  else if(any(grepl("muscle", x))) { y <- "muscle" }
  else if(any(grepl("development", x))) { y <- "development" }
  else if(any(grepl("morphogenesis", x))) { y <- "morphogenesis" }
  else if(any(grepl("extracellular matrix organization", x))) { y <- "ECM organization" }
  else if(any(grepl("histone modification", x))) { y <- "histone modification" }
  else if(any(grepl("posttranscriptional regulation", x))) { y <- "posttranscriptional regulation" }
  else if(any(grepl("regulation of gene expression", x))) { y <- "gene expression regulation" }
  else if(any(grepl("RNA splicing", x))) { y <- "RNA splicing" }
  else if(any(grepl("translation", x))) { y <- "translation" }
  else if(any(grepl("protein folding", x))) { y <- "protein folding" }
  else if(any(grepl("protein targeting", x))) { y <- "protein targeting" }
  #else if(any(grepl("protein transport", x))) { y <- "protein transport" }
  else if(any(grepl("vesicle transport", x))) { y <- "vesicle transport" }
  else if(any(grepl("respiration", x))) { y <- "respiration" }
  else if(any(grepl("antigen processing and presentation", x))) { y <- "immune" }
  #else if(any(grepl("regulation of cell migration", x))) { y <- "regulation of cell migration" }
  else if(any(grepl("regulation of Wnt signaling pathway", x))) { y <- "regulation of Wnt signaling pathway" }
  else if(any(grepl("cell cycle", x))) { y <- "cell cycle" }
  else if(any(grepl("ribonucleoprotein complex", x))) { y <- "ribonucleoprotein complex" }
  else if(any(grepl("acetylation", x))) { y <- "acetylation" }
  else if(any(grepl("metabolic", x))) { y <- "metabolism" }
  else (y <- "other")
  return(y)
})
###
mdmap_ds["MD25"] <- "translation"
mdmap_ds["MD30"] <- "ECM organization"
mdmap_ds["MD82"] <- "ECM organization"
mdmap_ds["MD91"] <- "ECM organization"
mdmap_ds["MD100"] <- "development"
mdmap_ds[c("MD97","MD123","MD126","MD135","MD136","MD146")] <- "immune"
###
mdmap_ds <- data.frame(mdmap_ds = Hmisc::capitalize(mdmap_ds), stringsAsFactors = F)
mdmap_ds_tmp <- data.frame(row.names = setdiff(rownames(md_map), rownames(mdmap_ds)), mdmap_ds = rep("Unknown", length(setdiff(rownames(md_map), rownames(mdmap_ds)))), stringsAsFactors = F)
mdmap_ds <- rbind(mdmap_ds, mdmap_ds_tmp)
colnames(mdmap_ds) <- "class"
mdmap_ds <- mdmap_ds[rownames(md_map), , drop = F]
mdmap_ds$class <- factor(mdmap_ds$class, levels = c(setdiff(sort(unique(mdmap_ds$class)), c("Other", "Unknown")), c("Other", "Unknown")))

# color
col_fun <- circlize::colorRamp2(c(-0.3, 0, 0.3), c("blue", "white", "red"))
col_fun_maxIoU <- circlize::colorRamp2(c(0, 0.5), c("#F7FCF5", "#41AB5D"))
# value
mdmap_mat <- md_map
colnames(mdmap_mat) <- mdmap_ct
# left anno
size_mat <- md_map_size
# right anno
enrich_mat <- as.matrix(mes$merged$cl_list_ftd[, c("enrich_TF", "enrich_PPI", "enrich_GO", "enrich_KEGG")])
colnames(enrich_mat) <- gsub("^enrich_", "", colnames(enrich_mat))
# top anno
load("03-expression/merged/cellCluster/ct_color.RData")
ts_DF <- data.frame(color = ct_color, stringsAsFactors = F)
ts_color <- ts_DF[mdmap_ts, ]
names(ts_color) <- mdmap_ts
# maxIoU (MsigDB)
maxIoU_mat <- read.table(file = "03-expression/merged/mdIoU/geneModule_maxIoU.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
colnames(maxIoU_mat) <- "Max.IoU"
# des
des_mat <- mdmap_ds
# expr
ex_mat <- mes$merged$avgExpr

## sort by avgCor
od <- order(rowMeans(mdmap_mat), decreasing = T)
mdmap_mat <- mdmap_mat[od, ]
size_mat <- size_mat[od, ]
enrich_mat <- enrich_mat[od, ]
maxIoU_mat <- maxIoU_mat[od, ]
des_mat <- des_mat[od, , drop = F]
ex_mat <- ex_mat[od, ]
##

# sort tissue
cellTypeMeta <- read.table(file = "cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
ts_ordered <- intersect(Hmisc::capitalize(unique(cellTypeMeta$tissue)), mdmap_ts)
#

# plot
pdf(paste0(OUT, "/module_map.pdf"), width = 8, height = 16)

mdid_case_ids <- c("MD51", "MD117", "MD91")
ht1 <- Heatmap(matrix = mdmap_mat, col = col_fun, name = "Correlation", 
               row_title_side = "right", row_title_rot = 0, row_title_gp = gpar(fontsize = 12), 
               column_title = "Cell type", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
               cluster_rows = F, row_split = des_mat$class, row_gap = unit(0.3, "mm"), border = T, 
               show_row_names = F, row_names_gp = gpar(fontsize = 3), 
               show_column_names = T, column_names_gp = gpar(fontsize = 12), column_names_rot = 45, 
               top_annotation = columnAnnotation(Organ = mdmap_ts, col = list(Organ = ts_color), simple_anno_size = unit(0.5, "cm"), 
                                                 annotation_name_gp = gpar(fontsize = 12), show_legend = F), 
               left_annotation = rowAnnotation(Log10.Size = anno_barplot(log10(size_mat$size), width = unit(2.5, "cm"), border = F, gp = gpar(fill = "#a1ccf7", col = NA), axis_param = list(gp = gpar(fontsize = 12))), 
                                               foo = anno_mark(at = match(mdid_case_ids, rownames(mdmap_mat)), labels = mdid_case_ids, side = "left", labels_gp = gpar(fontsize = 12), link_width = unit(3, "mm")), 
                                               annotation_name_gp = gpar(fontsize = 12), annotation_name_offset = unit(0.6, "cm")), 
               right_annotation = rowAnnotation(Enrichment = enrich_mat, Max.IoU = maxIoU_mat, col = list(Enrichment = c("TRUE" = "skyblue", "FALSE" = "grey90"), Max.IoU = col_fun_maxIoU), 
                                                annotation_name_gp = gpar(fontsize = 12), show_legend = F), 
               show_heatmap_legend = F)

#ht2 <- Heatmap(matrix = log10(ex_mat + 1))
lgd0 <- Legend(col_fun = col_fun, title = "Correlation", title_gp = gpar(fontsize = 12), title_gap = unit(2, "mm"), labels_gp = gpar(fontsize = 12), legend_height = unit(2, units = "cm"))
lgd1 <- Legend(at = ts_ordered, title = "Organ", legend_gp = gpar(fill = ct_color[ts_ordered]), title_gp = gpar(fontsize = 12), title_gap = unit(2, "mm"), labels_gp = gpar(fontsize = 12))
lgd2 <- Legend(at = c("True", "False"), title = "Enrichment", legend_gp = gpar(fill = c("skyblue", "grey90")), title_gp = gpar(fontsize = 12), title_gap = unit(2, "mm"), labels_gp = gpar(fontsize = 12))
lgd3 <- Legend(col_fun = col_fun_maxIoU, title = "Maximum IoU", title_gp = gpar(fontsize = 12), title_gap = unit(2, "mm"), labels_gp = gpar(fontsize = 12), legend_height = unit(2, units = "cm"))
lgd <- packLegend(lgd0, lgd1, lgd2, lgd3, direction = "horizontal", column_gap = unit(0.725, "cm"))
draw(ht1, row_title = "Gene module", row_title_gp = gpar(fontsize = 14), padding = unit(c(95, 5.5, 5.5, 5.5), units = "points"))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.01, "npc"), just = c("center", "bottom"))

dev.off()

## module case
# 0) global plot (specific cell type)
pdf(file = paste0(OUT, "/module_global.pdf"), width = 5.5, height = 5.5, useDingbats = F)

do_plotCorHeatmap(res_in = res, ctype1 = "Stomach.Fibro-ADAM28", ctype2 = NULL, mdid = unique(res$`Stomach.Fibro-ADAM28`$cl_list_ftd$cluster), show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + 
  theme(legend.position = "bottom", legend.justification = "center") + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(4, units = "cm")))

dev.off()

cor_ATAC <- read.table(file = paste0(OUT, "/Stomach.Fibro-ADAM28_ga_Jaccard.txt"), header = T, sep = "\t", stringsAsFactors = F)
dim(cor_ATAC)
cor_ATAC <- as.matrix(cor_ATAC)
res_ATAC <- res["Stomach.Fibro-ADAM28"]
res_ATAC$`Stomach.Fibro-ADAM28`$cor_cld <- cor_ATAC

pdf(file = paste0(OUT, "/module_global_ATAC.pdf"), width = 5.5, height = 5.5, useDingbats = F)

do_plotCorHeatmap(res_in = res_ATAC, ctype1 = "Stomach.Fibro-ADAM28", ctype2 = NULL, mdid = unique(res$`Stomach.Fibro-ADAM28`$cl_list_ftd$cluster), show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + 
  theme(legend.position = "bottom", legend.justification = "center") + 
  labs(fill = "Jaccard index") + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(4, units = "cm")))

dev.off()

# *) experiment validation
# MD91
pdf(file = paste0(OUT, "/module_case_MD91.pdf"), width = 8, height = 4, useDingbats = F)

ctype1_case <- "Small intestine.Fibro-COL6A5"
ctype2_case <- "Pancreas.Fibro-PAMR1+SOX6+"
mdid_case <- "MD91"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
  theme(legend.box.margin = margin(l = -12)) + ylab("Module genes")
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 0.9, do.print = F, do.return = T) + xlab("GO terms\n")
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
#do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)

dev.off()

# MD203
pdf(file = paste0(OUT, "/module_case_MD203.pdf"), width = 8, height = 4, useDingbats = F)

ctype1_case <- "Small intestine.SM-Visceral"
ctype2_case <- "Stomach.Fibro-IRF1"
mdid_case <- "MD203"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
  theme(legend.box.margin = margin(l = -12)) + ylab("CGM Genes")
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 0.8, do.print = F, do.return = T) + xlab("GO terms\n")
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
#do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)

dev.off()

# 1) same tissue, different cell types (muscle)
pdf(file = paste0(OUT, "/module_case_MD117.pdf"), width = 8, height = 4, useDingbats = F)

ctype1_case <- "Small intestine.SM-Visceral"
ctype2_case <- "Small intestine.Fibro-COL6A5"
mdid_case <- "MD117"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
  theme(legend.box.margin = margin(l = -12)) + ylab("Module genes")
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 0.5, do.print = F, do.return = T) + xlab("GO terms\n") + 
  scale_y_continuous(breaks = c(0, 5, 10))
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
#do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)

dev.off()

# KEGG
# g_cand <- strsplit(subset(mes$merged$enrich_KEGG, cluster == "MD117" & ID == "hsa04510", "geneID", drop = T), "/")[[1]]
# gd <- rep(-1, length(g_cand))
# names(gd) <- g_cand
# pathview(gene.data = gd, gene.idtype = "SYMBOL", pathway.id = "hsa04510", species = "hsa", out.suffix = "XXX")

# PPI
# https://string-db.org/cgi/network.pl?taskId=vO1dZKWPO4ZX

# 2) different tissues (development)
pdf(file = paste0(OUT, "/module_case_MD45.pdf"), width = 8, height = 4, useDingbats = F)

ctype1_case <- "Ovary.Granulosa-R-Al"
ctype2_case <- "Small intestine.T"
mdid_case <- "MD45"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
  theme(legend.box.margin = margin(l = -12))
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
#do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)

dev.off()

# 3) ubiquent
# MD51
pdf(file = paste0(OUT, "/module_case_MD51.pdf"), width = 8, height = 4, useDingbats = F)

ctype1_case <- "Stomach.Fibro-FBLN1"
ctype2_case <- "Pancreas.Fibro-PAMR1+SOX6+"
mdid_case <- "MD51"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
  theme(legend.box.margin = margin(l = -12)) + ylab("Module genes")
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 0.9, do.print = F, do.return = T) + xlab("GO terms\n") + 
  scale_y_continuous(breaks = c(0, 5, 10))
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
#do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)

dev.off()

# KEGG
# g_cand <- strsplit(subset(mes$merged$enrich_KEGG, cluster == "MD51" & ID == "hsa04010", "geneID", drop = T), "/")[[1]]
# gd <- rep(-1, length(g_cand))
# names(gd) <- g_cand
# pathview(gene.data = gd, gene.idtype = "SYMBOL", pathway.id = "hsa04010", species = "hsa", out.suffix = "XXX")

# PPI
# https://string-db.org/cgi/network.pl?taskId=uGoHw0StpUZs

# MD16
pdf(file = paste0(OUT, "/module_case_MD16.pdf"), width = 8, height = 4, useDingbats = F)

ctype1_case <- "Stomach.Fibro-FBLN1"
ctype2_case <- "Pancreas.Fibro-PAMR1+SOX6+"
mdid_case <- "MD16"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
  theme(legend.box.margin = margin(l = -12))
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
#do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)

dev.off()

# DEGs
DEGs <- read.table(file = "03-expression/merged/cellCluster/expr.markers_ftd_Sm.Fibro-COL6A5_Pa.Fibro-PAMR1+SOX6+.txt", header = T, sep = "\t", stringsAsFactors = F)
DEGs <- DEGs[order(DEGs$isMD91), ]

pdf(file = paste0(OUT, "/module_case_DEGs.pdf"), width = 5.5, height = 4, useDingbats = F)

ctype1_case <- "Small intestine.Fibro-COL6A5"
ctype2_case <- "Pancreas.Fibro-PAMR1+SOX6+"
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mgid = DEGs$gene, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(ctype1_case)

dev.off()

# pdf(file = paste0(OUT, "/module_case_MD76.pdf"), width = 8, height = 4)
# 
# ctype1_case <- "Stomach.Fibro-FBLN1"
# ctype2_case <- "Small intestine.T"
# mdid_case <- "MD76"
# do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
# do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
#   theme(legend.box.margin = margin(l = -12))
# do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 3)
# do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
# do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
# #do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)
# 
# dev.off()

# pdf(file = paste0(OUT, "/module_case_MD101.pdf"), width = 8, height = 4)
# 
# ctype1_case <- "Stomach.Fibro-FBLN1"
# ctype2_case <- "Small intestine.T"
# mdid_case <- "MD101"
# do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T) + ylab(paste0(mdid_case, "\n\n", ctype1_case))
# do_plotExprHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, mask0 = F, show_rownames = F, fontsize_row = 10, short_name = T, asp.ratio = 2.75, do.print = F, do.return = T) + 
#   theme(legend.box.margin = margin(l = -12))
# do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "GO", main = NULL, asp.ratio = 3)
# do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "KEGG", main = NULL, asp.ratio = 3)
# do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "TF", main = NULL, bar.col = "aquamarine2", asp.ratio = 3)
# #do_plotEnrich(res_in = mes, ctype = "merged", mdid = mdid_case, db = "miRNA", main = NULL)
# 
# dev.off()

# 4) just for showing
pdf(file = paste0(OUT, "/module_case_simple.pdf"), width = 6, height = 6)

ctype1_case <- "Small intestine.SM-Visceral"
ctype2_case <- "Small intestine.Fibro-COL6A5"
mdid_case <- c("MD51", "MD117")
do_plotCorHeatmap(res_in = res, ctype1 = ctype1_case, ctype2 = ctype2_case, mp_in = cl_rmdup_LS, mpid = mdid_case, show_rownames = F, fontsize_row = 10, do.print = F, do.return = T, show.legend = F, do_highlights = T) + 
  xlab(NULL) + ylab(NULL)

dev.off()

# expr ~ module
hc_module <- hclust(dist(t(md_map)), method = "complete")
plot(hc_module)

expr_map <- sapply(res, function(x) { y <- rowMeans(x$data); return(y) })
hc_expr <- hclust(dist(t(expr_map)), method = "complete")
plot(hc_expr)

dend12 <- dendlist(as.dendrogram(hc_expr), as.dendrogram(hc_module))
dl12 <- dend12 %>% untangle(method = "step2side")

pdf(file = paste0(OUT, "/tanglegram_expr_module.pdf"), width = 6, height = 3.5)
tanglegram(dl12, lab.cex = 1.2, cex_main = 0.8, cex_main_left = 1.25, cex_main_right = 1.25, axes = F, 
           lwd = 1, edge.lwd = 1, common_subtrees_color_branches = TRUE, common_subtrees_color_lines = TRUE,
           columns_width = c(18, 2, 18), 
           main_left = "Expression", main_right = "Module",
           margin_inner= 16, margin_bottom = 0, 
           main = paste("Entanglement =", round(entanglement(dl12), 3)))
dev.off()

# gene type composition ----
cl_list_ftd_alt <- mes$merged$cl_list_ftd[, c("cluster", "protein_coding", "lncRNA", "sncRNA", "pseudogene", "Ig/TCR")]
### re-class
cl_list_ftd_alt$protein_coding <- cl_list_ftd_alt$protein_coding + cl_list_ftd_alt$`Ig/TCR`
cl_list_ftd_alt$lncRNA <- cl_list_ftd_alt$lncRNA + cl_list_ftd_alt$pseudogene
###
cl_list_ftd_gstat_melted <- reshape2::melt(cl_list_ftd_alt[, 1:4], id.vars = "cluster")
###
#md_ordered <- rownames(mes$merged$avgCor)[order(apply(mes$merged$avgCor, 1, max))]
md_ordered <- cl_list_ftd_alt$cluster[order(cl_list_ftd_alt$protein_coding, decreasing = T)]
###
cl_list_ftd_gstat_melted$cluster <- factor(cl_list_ftd_gstat_melted$cluster, levels = md_ordered)
cl_list_ftd_gstat_melted$variable <- factor(cl_list_ftd_gstat_melted$variable, levels = c(levels(cl_list_ftd_gstat_melted$variable), "TF"))
levels(cl_list_ftd_gstat_melted$variable) <- c("PCG", "LncRNA", "SncRNA                          ", "TF")

pdf(paste0(OUT, "/geneTypeComposition.pdf"), width = 6.5, height = 4)

ggplot(cl_list_ftd_gstat_melted, aes(x = cluster, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = T)) + 
  geom_bar(data = mes$merged$cl_list_ftd, aes(y = tf), fill = "brown1", stat = "identity") + 
  #coord_flip(clip = "off") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.direction = "horizontal", legend.position = "bottom", legend.justification = "center", legend.box.margin = margin(l = -5, t = -16, b = -6)) + 
  theme(plot.margin = margin(10,7,7,7)) + 
  labs(fill = NULL) + ylab("Gene type composition") + 
  guides(fill = guide_legend(nrow = 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  scale_fill_manual(values = c("hotpink", "dodgerblue", "lightskyblue", "grey60", "grey90", "brown1")[c(1:3,6)], drop = F)

dev.off()

# genomic distances ----
write.table(x = unique(mes$merged$cl_table_ftd$gene), file = paste0(OUT, "/genes_in_module.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
system(paste0("perl do_geneDist.pl ", OUT, "/genes_in_module.txt > ", OUT, "/geneDist.txt"))
geneDist <- read.table(paste0(OUT, "/geneDist.txt"), header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
dim(geneDist)
colnames(geneDist) <- c("gene1", "gene2", "ds")
geneDist_MT <- acast(geneDist, gene1 ~ gene2, value.var = "ds")
geneDist_stat <- t(sapply(cl_rmdup_LS, function(x) {
  xds <- geneDist_MT[x, x]
  xds <- xds[upper.tri(xds)]
  ds_rat <- sum(! is.na(xds)) / length(xds)
  ds_avg <- mean(xds, na.rm = T)
  if(is.na(ds_avg)) { ds_avg <- Inf }
  y <- c(ds_rat, ds_avg)
  return(y)
}))
colnames(geneDist_stat) <- c("ds_rat", "ds_avg")
geneDist_stat <- as.data.frame(geneDist_stat)
geneModule_stat <- merge(mes$merged$cl_list_ftd, geneDist_stat, by.x = "cluster", by.y = 0, sort = F)

pdf(paste0(OUT, "/geneDist.pdf"), width = 6.5, height = 4)

ggplot(geneModule_stat, aes(x = ds_rat, y = ds_avg / 1e6)) + 
  geom_point(aes(size = size, color = pmin(apply(mes$merged$avgCor, 1, max), 0.3)), alpha = 0.6) + 
  ggrepel::geom_text_repel(data = subset(geneModule_stat, ds_rat > 0.5 | is.infinite(ds_avg)), aes(label = cluster), box.padding = 0.5, point.padding = 0.6) + 
  coord_cartesian(clip = "off") + xlab("Ratio of gene pairs in the same chromosome") + ylab("Average genomic distance (Mb)") + 
  labs(color = "Max. correlation", size = "Size") + guides(size = guide_legend(override.aes = list(color = "red"), order = 1)) + 
  scale_color_gradientn(colors = colorRampPalette(c("blue", "white", "red"))(100), limits = c(-0.3, 0.3), breaks = c(-0.3, -0.3 / 2, 0, 0.3 / 2, 0.3))

dev.off()

## avgCor ~ TF num
md_avgExpr <- mes$merged$avgExpr
md_avgCor <- mes$merged$avgCor
md_avgCor_max <- as.data.frame(apply(md_avgCor, 1, max))
colnames(md_avgCor_max) <- "avgCor_max"
md_TFnum <- as.data.frame(table(mes$merged$enrich_TF$cluster))
colnames(md_TFnum) <- c("module_id", "TF_num")
md_TFnum$module_id <- as.character(md_TFnum$module_id)
md_DF <- merge(md_avgCor_max, md_TFnum, by.x = 0, by.y = "module_id", sort = F, all.x = T)
md_DF[is.na(md_DF$TF_num), "TF_num"] <- 0
ggplot(md_DF, aes(x = avgCor_max, y = TF_num)) + geom_point()
##

# avgCor ~ enrich TF
md_enrichTFbyAvgCor_LS <- split(mes$merged$cl_list_ftd$enrich_TF, cut(apply(mes$merged$avgCor, 1, max), breaks = c(0, 0.1, 0.2, 0.3, 1)))
md_enrichTFbyAvgCor_MT <- sapply(md_enrichTFbyAvgCor_LS, function(x) { c(sum(! x), sum(x), sum(x) / length(x)) })
rownames(md_enrichTFbyAvgCor_MT) <- c("without_TF", "with_TF", "ratio")
fisher.test(md_enrichTFbyAvgCor_MT[1:2, c(1, 4)], alternative = "greater")$p.value
md_enrichTFbyAvgCor <- data.frame(avgCor_max = names(md_enrichTFbyAvgCor_LS), ratio = md_enrichTFbyAvgCor_MT["ratio", ])

pdf(paste0(OUT, "/md_enrichTFbyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichTFbyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched TFs") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

# avgCor ~ enrich PPI
md_enrichPPIbyAvgCor_LS <- split(mes$merged$cl_list_ftd$enrich_PPI, cut(apply(mes$merged$avgCor, 1, max), breaks = c(0, 0.1, 0.2, 0.3, 1)))
md_enrichPPIbyAvgCor_MT <- sapply(md_enrichPPIbyAvgCor_LS, function(x) { c(sum(! x), sum(x), sum(x) / length(x)) })
rownames(md_enrichPPIbyAvgCor_MT) <- c("without_PPI", "with_PPI", "ratio")
fisher.test(md_enrichPPIbyAvgCor_MT[1:2, c(1, 4)], alternative = "greater")$p.value
md_enrichPPIbyAvgCor <- data.frame(avgCor_max = names(md_enrichPPIbyAvgCor_LS), ratio = md_enrichPPIbyAvgCor_MT["ratio", ])

pdf(paste0(OUT, "/md_enrichPPIbyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichPPIbyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched PPI") + #ggtitle("STRING") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

pdf(paste0(OUT, "/md_enrichSTRINGbyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichPPIbyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched PPI") + ggtitle("STRING") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

# avgCor ~ enrich GO
md_enrichGObyAvgCor_LS <- split(mes$merged$cl_list_ftd$enrich_GO, cut(apply(mes$merged$avgCor, 1, max), breaks = c(0, 0.1, 0.2, 0.3, 1)))
md_enrichGObyAvgCor_MT <- sapply(md_enrichGObyAvgCor_LS, function(x) { c(sum(! x), sum(x), sum(x) / length(x)) })
rownames(md_enrichGObyAvgCor_MT) <- c("without_GO", "with_GO", "ratio")
fisher.test(md_enrichGObyAvgCor_MT[1:2, c(1, 4)], alternative = "greater")$p.value
md_enrichGObyAvgCor <- data.frame(avgCor_max = names(md_enrichGObyAvgCor_LS), ratio = md_enrichGObyAvgCor_MT["ratio", ])

pdf(paste0(OUT, "/md_enrichGObyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichGObyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched GO") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

# avgCor ~ enrich KEGG
md_enrichKEGGbyAvgCor_LS <- split(mes$merged$cl_list_ftd$enrich_KEGG, cut(apply(mes$merged$avgCor, 1, max), breaks = c(0, 0.1, 0.2, 0.3, 1)))
md_enrichKEGGbyAvgCor_MT <- sapply(md_enrichKEGGbyAvgCor_LS, function(x) { c(sum(! x), sum(x), sum(x) / length(x)) })
rownames(md_enrichKEGGbyAvgCor_MT) <- c("without_KEGG", "with_KEGG", "ratio")
fisher.test(md_enrichKEGGbyAvgCor_MT[1:2, c(1, 4)], alternative = "greater")$p.value
md_enrichKEGGbyAvgCor <- data.frame(avgCor_max = names(md_enrichKEGGbyAvgCor_LS), ratio = md_enrichKEGGbyAvgCor_MT["ratio", ])

pdf(paste0(OUT, "/md_enrichKEGGbyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichKEGGbyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched KEGG") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

# avgCor ~ enrich HIPPIE (20200804)
md_enrichHIPPIEbyAvgCor_LS <- split(mes$merged$cl_list_ftd$enrich_HIPPIE, cut(apply(mes$merged$avgCor, 1, max), breaks = c(0, 0.1, 0.2, 0.3, 1)))
md_enrichHIPPIEbyAvgCor_MT <- sapply(md_enrichHIPPIEbyAvgCor_LS, function(x) { c(sum(! x), sum(x), sum(x) / length(x)) })
rownames(md_enrichHIPPIEbyAvgCor_MT) <- c("without_HIPPIE", "with_HIPPIE", "ratio")
fisher.test(md_enrichHIPPIEbyAvgCor_MT[1:2, c(1, 4)], alternative = "greater")$p.value
md_enrichHIPPIEbyAvgCor <- data.frame(avgCor_max = names(md_enrichHIPPIEbyAvgCor_LS), ratio = md_enrichHIPPIEbyAvgCor_MT["ratio", ])

pdf(paste0(OUT, "/md_enrichHIPPIEbyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichHIPPIEbyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "**", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched PPI") + ggtitle("HIPPIE") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

# avgCor ~ enrich HuRI (20200804)
md_enrichHuRIbyAvgCor_LS <- split(mes$merged$cl_list_ftd$enrich_HuRI, cut(apply(mes$merged$avgCor, 1, max), breaks = c(0, 0.1, 0.2, 0.3, 1)))
md_enrichHuRIbyAvgCor_MT <- sapply(md_enrichHuRIbyAvgCor_LS, function(x) { c(sum(! x), sum(x), sum(x) / length(x)) })
rownames(md_enrichHuRIbyAvgCor_MT) <- c("without_HuRI", "with_HuRI", "ratio")
fisher.test(md_enrichHuRIbyAvgCor_MT[1:2, c(1, 4)], alternative = "greater")$p.value
md_enrichHuRIbyAvgCor <- data.frame(avgCor_max = names(md_enrichHuRIbyAvgCor_LS), ratio = md_enrichHuRIbyAvgCor_MT["ratio", ])

pdf(paste0(OUT, "/md_enrichHuRIbyAvgCor.pdf"), width = 3.25, height = 4.5)

ymax <- 1 * 1.2
sigBar <- data.frame(avgCor_max = rep(c("(0,0.1]", "(0.3,1]"), each = 2), ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(md_enrichHuRIbyAvgCor, aes(x = avgCor_max, y = ratio)) + geom_bar(stat = "identity", fill = "cornflowerblue") + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 2.5, y = ymax * (1 - 0.05), label = "*", col="red", size = 5.5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  xlab("Max. correlation") + ylab("Ratio of modules with enriched PPI") + ggtitle("HuRI") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10,7,4,7))

dev.off()

# 6. gene module metatable ----
md_meta <- data.frame(size_mat, enrich_mat, des_mat, row.names = NULL, stringsAsFactors = F, check.names = F)

mdmap_mat_new <- mdmap_mat
colnames(mdmap_mat_new) <- colnames(mes$merged$avgCor)
colnames(mdmap_mat_new) <- paste0("avgCor_", colnames(mdmap_mat_new))

ex_mat_new <- ex_mat
colnames(ex_mat_new) <- colnames(mes$merged$avgExpr)
colnames(ex_mat_new) <- paste0("avgExpr_", colnames(ex_mat_new))

md_meta <- cbind(md_meta, mdmap_mat_new, ex_mat_new)
# sort by des
md_meta <- do.call("rbind", split(md_meta, md_meta$class))
rownames(md_meta) <- NULL

write.table(x = md_meta, file = paste0(OUT, "/geneModule_metatable.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# cell type order
ct_ordered <- gsub("^avgCor_", "", colnames(mdmap_mat_new)[c(7,2,1,9,3,6,4,5,10,8)])
write.table(x = ct_ordered, file = paste0(OUT, "/cellType_ordered.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# module to gene
md_genes <- mes$merged$cl_table_ftd[, 1:2]
write.table(x = md_genes, file = paste0(OUT, "/geneModule_genes.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# cor mat
res_cor_inMD <- lapply(res, function(x) {
  g_sub <- intersect(rownames(x$cor_cld), mes$merged$cl_table_ftd$gene)
  x$cor_cld <- x$cor_cld[g_sub, g_sub]
  x$cor_cld <- expandCorMat(cand = x$cor_cld, ref = matrix(NA, 
                                                           nrow = length(unique(mes$merged$cl_table_ftd$gene)), 
                                                           ncol = length(unique(mes$merged$cl_table_ftd$gene)), 
                                                           dimnames = list(unique(mes$merged$cl_table_ftd$gene), unique(mes$merged$cl_table_ftd$gene))))
  y <- x[["cor_cld"]]
  return(y)
})
saveRDS(res_cor_inMD, file = paste0(OUT, "/res_cor_inMD.rds"))

# expr mat
res_expr_inMD <- lapply(res, function(x) {
  set.seed(1)
  x$data <- x$data[rownames(x$data) %in% mes$merged$cl_table_ftd$gene, sample(ncol(x$data), 500)]
  x$data <- x$data[, order(colMeans(sweep(x$data, 1, rowSums(x$data), "/")))]
  y <- x[["data"]]
  return(y)
})
saveRDS(res_expr_inMD, file = paste0(OUT, "/res_expr_inMD.rds"))

# 7. if TF control the holdOn of module ----
# ID mapping
id_mapping_LS <- list(NCOR="NCOR1", PU="SPI1", EKLF="KLF1", JARID1A="KDM5A", 
                      OCT4="POU5F1", EST1="ETS1", GABP="GABPA", TCFCP2L1="TFCP2L1", 
                      TCFAP2C="TFAP2C", RXR="RXRA", P300="EP300", STAT5="STAT5A", 
                      SFPI1="SPI1", DPY="DPY30")
id_mapping_DF <- melt(id_mapping_LS)[, c(2,1)]
colnames(id_mapping_DF) <- c("old", "new")
id_mapping_DF$old <- as.character(id_mapping_DF$old)
id_mapping_DF$new <- as.character(id_mapping_DF$new)

calckeyTF <- function(cmpMd_DF, enrichment_ftd_LS, ctype1, ctype2, diff_DF, corm_DF, do_plot1 = F, do_plot2 = F, fontsize_row) {
  er1 <- enrichment_ftd_LS[[ctype1]]
  et1 <- subset(er1, database == "ChEA_2016" | database == "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  et1$TF <- gsub("_.*", "", et1$Term)
  ets1 <- merge(subset(cmpMd_DF, holdOn ==0), et1, by.x = 0, by.y = "cluster", sort = F)
  colnames(ets1)[1] <- "cluster"
  # diff gene
  diffGene_DF <- diffGene(expr_data, cellMetaData, ctype1, ctype2)
  ets1$isDftb <- ets1$TF %in% rownames(diffGene_DF)
  # check missed ID
  geneID <- read.table("~/lustre/06-Human_cell_atlas/Data/gene_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)[, 2]
  ets1$isAnno <- ets1$TF %in% geneID
  etsm1 <- merge(ets1, id_mapping_DF, by.x = "TF", by.y = "old", all.x = T, sort = F)
  
  cat("Unknown genes:", nrow(unique(subset(etsm1, isDftb==F, "TF"))), "\n")
  cat("  No expression:", nrow(unique(subset(etsm1, isDftb==F & isAnno==T, "TF"))), "\n")
  cat("  Not in gene list:", nrow(unique(subset(etsm1, isDftb==F & isAnno==F, "TF"))), "\n")
  cat("  Rescued:", nrow(unique(subset(etsm1, isDftb==F & isAnno==F & (! is.na(new)), "TF"))), "\n")
  
  etsm1[etsm1$isAnno, "new"] <- etsm1[etsm1$isAnno, "TF"]
  etsm1$fnAnno <- etsm1$new %in% geneID
  etsd1 <- merge(etsm1, diffGene_DF, by.x = "new", by.y = 0, sort = F)
  etsd1 <- subset(etsd1, logFC<0 & qvalue<0.05)
  
  cat("------------\n")
  cat("Total module:", nrow(cmpMd_DF), "\n")
  cat("  Show high correlation:", sum(cmpMd_DF$holdOn==1), "\n")
  cat("  Show low correlation:", sum(cmpMd_DF$holdOn==0), "\n")
  cat("    Contain at least 1 common TF:", nrow(unique(subset(etsm1, holdOn==0, "cluster"))), "\n")
  cat("      Contain at least 1 decreased common TF:", nrow(unique(subset(etsd1, holdOn==0, "cluster"))), "\n")
  
  if(do_plot1) {
    pdf(paste0("03-expression/merged/geneModule/", samplingPos, "/calckeyTF_", ctype1, "_", ctype2, ".pdf"), height = 5.5, width = 3)
    anno_row_DF <- data.frame(row.names = rownames(diff_DF), 
                              Decline = abs(diff_DF[, "holdOn"]-1), 
                              TFsupport = as.numeric( (rownames(corm_DF)%in%etsd1$cluster) & (diff_DF[, "holdOn"] == 0) ))
    anno_row_DF$Decline <- factor(anno_row_DF$Decline)
    anno_row_DF$TFsupport <- factor(anno_row_DF$TFsupport)
    p <- pheatmap(corm_DF, cluster_rows = F, cluster_cols = F, 
                  annotation_row = anno_row_DF[, 2:1], 
                  annotation_colors = list(Decline = c(`0` = "grey", `1` = "green"), TFsupport = c(`0` = "grey", `1` = "skyblue")), 
                  annotation_names_row = F, angle_col = 45, fontsize_row = fontsize_row, silent = T)
    print(p)
    dev.off()
  }
  return(etsd1)
}

do_diffMd <- function (cluster_table_LS, cor_LS, ctype1, ctype2) {
  ct1 <- cluster_table_LS[[ctype1]]
  ct2 <- cluster_table_LS[[ctype2]]
  cor1 <- cor_LS[[ctype1]][ct1$gene, ct1$gene]
  cor2 <- expandCorMat(cor_LS[[ctype2]], cor1)[ct1$gene, ct1$gene]  # use the genes in ct1
  # diff module
  diff_DF <- diffMd(cor1, cor2, ct1, cutoff = 0.3)
  
  # list the average pair-wise correlation for each module
  ctc1_ft <- factor(ct1$cluster, levels = unique(ct1$cluster))
  corm1 <- corMatStat(mt = cor1, cond = list(ctc1_ft))
  corm2 <- corMatStat(mt = cor2, cond = list(ctc1_ft))
  corm_DF <- data.frame(diag(corm1), diag(corm2))
  colnames(corm_DF) <- c(ctype1, ctype2)
  
  out <- list(diff_DF = diff_DF, corm_DF = corm_DF)
  return(out)
}

diffMd_res <- do_diffMd(cluster_table_LS, cor_LS, "COL17A1+", "collagen+")

calckeyTF_res <- calckeyTF(cmpMd_DF1, enrichment_ftd_LS, "COL17A1+", "collagen+", diffMd_res$diff_DF, diffMd_res$corm_DF, 
                           do_plot1 = T, do_plot2 = F, 3)

# plot case network
plot_caseNT <- function(calckeyTF_res, md1) {
  library("igraph")
  library("qgraph")
  calckeyTF_sub <- subset(calckeyTF_res, cluster == md1)
  calckeyTF_nt_LS <- apply(as.matrix(calckeyTF_sub[, c("new", "Genes")]), 1, function(x) { 
    y <- matrix(NA, nrow = 0, ncol = 2)
    for(i in strsplit(x[2], ";")[[1]]) { 
      y <- rbind(y, c(x[1], i))
    }
    return(y)
  })
  calckeyTF_nt_DF <- as.data.frame(do.call("rbind", calckeyTF_nt_LS), stringsAsFactors = F)
  g <- graph_from_data_frame(d = calckeyTF_nt_DF, directed = T)
  layout <- layout.reingold.tilford(g, circular=T)
  plot.igraph(g, layout = layout_with_fr)
}


# stat for all databases
#enrichStat <- function(cluster_table_LS, enrichment_ftd_LS, ctype1, ctype2, do_plot1 = F) {
eu_all <- list()
do_plot1 = T
pdf(paste0("03-expression/merged/geneModule/", samplingPos, "/enrichStat.pdf"), height = 4.5, width = 5.5)
for(ctype1 in names(cluster_table_LS)) {
  print(ctype1)
  ct1 <- cluster_table_LS[[ctype1]]
  er1 <- enrichment_ftd_LS[[ctype1]][, 1:9]
  eu1 <- unique(er1[, c("cluster", "database")])
  eu1_DF <- data.frame(table(eu1$database), unanno=length(unique(ct1$cluster))-as.numeric(table(eu1$database)))
  colnames(eu1_DF) <- c("annotation", "annotated", "unannotated")
  eu_all[[ctype1]] <- data.frame(ct = ctype1, eu1_DF[, 1:2], stringsAsFactors = F)
  if(do_plot1) {
    eu1_melted <- melt(eu1_DF, id.vars = "annotation", variable.name = "type", value.name = "number")
    ### change label
    levels(eu1_melted$annotation)
    levels(eu1_melted$annotation)[2] <- "ENCODE_ChEA_TF"
    levels(eu1_melted$annotation)
    ###
    eu1_melted$type <- factor(eu1_melted$type, levels = c("unannotated", "annotated"))
    
    #pdf(paste0("03-expression/merged/geneModule/", samplingPos, "/enrichStat_", ctype1, "_", ctype2, ".pdf"), height = 4.5, width = 5.5)
    p <- ggplot(eu1_melted, aes(x = annotation, y = number, fill = type)) + geom_bar(stat = "identity") + 
      scale_y_continuous(limits = c(0, length(unique(ct1$cluster)) * 1.05), expand = c(0, 0)) + 
      scale_fill_manual(values = c("grey", scales::hue_pal()(2)[2])) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + 
      ylab("Module number") + ggtitle(ctype1)
    print(p)
    #dev.off()
  }
}
dev.off()

names(eu_all) <- NULL
eu_all_DF <- do.call("rbind", eu_all)
md_stat_plus <- merge(md_stat, eu_all_DF, by.x = 0, by.y = "ct", sort = F)
colnames(md_stat_plus)[1] <- "ct"
write.table(x = md_stat_plus, file = paste0("03-expression/merged/geneModule/", samplingPos, "/md_stat_plus.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")

# X. save object ----
saveRDS(res, file = file.path(OUT, "res.rds"))
saveRDS(mes, file = file.path(OUT, "mes.rds"))
save.image(file = paste0(OUT, "/module.RData"))
