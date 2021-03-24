# check expr
setwd("~/lustre/06-Human_cell_atlas/pooled_data/KD/")

library("reshape2")
library("ggplot2")
library("cowplot")

samplingPos <- "."
OUT <- paste0("03-expression/merged/geneModule/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

# load expression
exprMatr_cellFiltered_CPM <- read.table(file = "03-expression/merged/filtering/UMIcount_cellFiltered_CPM.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F, comment.char = "")

# load meta
cellMetaData <- read.table(file = "../All/cell_metatable_filtered_plus.txt", header = T, sep = "\t", stringsAsFactors = F)
cellMetaData <- subset(cellMetaData, cell %in% colnames(exprMatr_cellFiltered_CPM), c("cell", "tissue", "samplingPos", "ident", "group"))
cellMetaData_tmp <- data.frame(cell = colnames(exprMatr_cellFiltered_CPM), tissue = "small intestine", samplingPos = NA, stringsAsFactors = F)
cellMetaData_tmp$ident <- gsub("^XC_", "", gsub("_[0-9][0-9]$", "", colnames(exprMatr_cellFiltered_CPM)))
cellMetaData_tmp$ident <- factor(cellMetaData_tmp$ident, levels = unique(cellMetaData_tmp$ident))
# change labels
levels(cellMetaData_tmp$ident)
levels(cellMetaData_tmp$ident)[1] <- "Control"
levels(cellMetaData_tmp$ident)
#
cellMetaData_tmp$group <- gsub("^XC_", "", gsub("-.*", "", cellMetaData_tmp$cell))
cellMetaData_tmp$group <- factor(cellMetaData_tmp$group, levels = unique(cellMetaData_tmp$group))
# changle labels
levels(cellMetaData_tmp$group)
levels(cellMetaData_tmp$group) <- c("Control", "siFOXL1", "siNR2F1")
levels(cellMetaData_tmp$group)
#
cellMetaData <- rbind(subset(cellMetaData, group != "X"), cellMetaData_tmp)
rm(cellMetaData_tmp)

# load cell type meta
cellTypeMeta <- read.table(file = "../All/cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
cellTypeMeta_sub <- subset(cellTypeMeta, tissue == cellMetaData$tissue[1])

# TF expression: barplot (by ident, cell level)
table(cellMetaData$ident)
res <- lapply(split(cellMetaData$cell, cellMetaData$ident), function(x) { y <- list(data = exprMatr_cellFiltered_CPM[, x]) })

tf_sub <- c("FOXL1", "NR2F1", "GAPDH")
expr_tf_LS <- lapply(names(res), function(x) {
  expr_tf_sub <- res[[x]]$data[tf_sub, ]
  set.seed(1)
  expr_tf_sub <- expr_tf_sub[, sample(1:ncol(expr_tf_sub), min(60, ncol(expr_tf_sub)))]
  expr_tf_sub$gene <- rownames(expr_tf_sub)
  expr_tf_sub_melted <- melt(expr_tf_sub, id.vars = "gene")
  expr_tf_sub_melted$cellType <- x
  return(expr_tf_sub_melted)
})
expr_tf_DF <- do.call("rbind", expr_tf_LS)
expr_tf_DF$gene <- factor(expr_tf_DF$gene, levels = tf_sub)
expr_tf_DF$cellType <- factor(expr_tf_DF$cellType, levels = unique(expr_tf_DF$cellType))
cd_col <- c("grey", RColorBrewer::brewer.pal(n = 6, name = "Dark2"))
names(cd_col) <- levels(expr_tf_DF$cellType)

pdf(file = paste0(OUT, "/geneExpr_barplot", ".pdf"), width = 8, height = 6, useDingbats = F)

ggplot(expr_tf_DF, aes(x = variable, y = log10(value + 1), fill = cellType)) + geom_bar(stat = "identity") + facet_grid(gene ~ .) + 
  #geom_vline(xintercept = c(70, 140), linetype = "dashed") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  xlab("Cell") + ylab(parse(text = "Log[10]~(CPM+1)")) + labs(fill = NULL) + 
  #ggtitle(paste("TFs in Small intestine")) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = cd_col)

dev.off()

# 2. comparison (by ident) ----
do_checkExpr <- function(TF_id) {
  cell_sub <- subset(cellMetaData, group %in% c("Control", paste0("si", TF_id)), "cell", drop = T)
  expr_sub <- as.data.frame(t(exprMatr_cellFiltered_CPM[TF_id, cell_sub, drop = F]))
  expr_sub$ident <- cellMetaData$ident[match(rownames(expr_sub), cellMetaData$cell)]
  colnames(expr_sub)[1] <- "CPM"
  stat_DF <- aggregate(CPM ~ ident, data = expr_sub, FUN = mean)
  stat_DF$fc <- stat_DF[, "CPM"] / stat_DF[1, "CPM"]
  stat_DF$pvalue <- aggregate(CPM ~ ident, data = expr_sub, FUN = function(x) { y <- wilcox.test(subset(expr_sub, ident == "Control", "CPM", drop = T), x, alternative = "greater")$p.value })[, 2]
  stat_DF$dtRatio <- aggregate(CPM ~ ident, data = expr_sub, FUN = function(x) { y <- sum(x > 0) / length(x) })[, 2]
  stat_DF$CPM_clm <- aggregate(CPM ~ ident, data = expr_sub, FUN = function(x) { set.seed(1); mean_cl_boot(x)[, 2] })[, 2]
  stat_DF$CPM_clp <- aggregate(CPM ~ ident, data = expr_sub, FUN = function(x) { set.seed(1); mean_cl_boot(x)[, 3] })[, 2]
  print(stat_DF)
  
  # plot
  ymax <- max(stat_DF$CPM_clp) * 1.2
  gp <- ggplot(stat_DF, aes(x = ident, y = CPM, fill = ident)) + 
    geom_col(show.legend = F) +
    geom_errorbar(aes(ymin = CPM_clm, ymax = CPM_clp), width = 0.2) + 
    geom_path(data = data.frame(ident = rep(unique(expr_sub$ident)[c(1, length(unique(expr_sub$ident)))], each = 2), CPM = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F), group = 1) + 
    annotate(geom = "text", x = (1 + length(unique(expr_sub$ident))) / 2, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
    xlab(NULL) + ylab("CPM") + ggtitle(TF_id) + 
    scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(0, ymax)) + 
    scale_fill_manual(values = cd_col[unique(expr_sub$ident)]) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
    theme(aspect.ratio = 1, plot.title = element_text(face = "plain")) + theme(aspect.ratio = 1.5)
  return(gp)
}

gp1 <- do_checkExpr(TF_id = "FOXL1")
gp2 <- do_checkExpr(TF_id = "NR2F1")

pdf(file = paste0(OUT, "/TF_expr_inKD_by_ident.pdf"), width = 5.5, height = 4)

gp1
gp2
plot_grid(gp1, gp2, align = "hv")

dev.off()

# 3. comparison (by group) ----
TF_id <- "FOXL1"
cell_sub <- subset(cellMetaData, group %in% c("Control", paste0("si", TF_id)), "cell", drop = T)
expr_sub <- as.data.frame(t(exprMatr_cellFiltered_CPM[TF_id, cell_sub, drop = F]))
expr_sub$group <- cellMetaData$group[match(rownames(expr_sub), cellMetaData$cell)]
colnames(expr_sub)[1] <- "CPM"
stat_DF <- aggregate(CPM ~ group, data = expr_sub, FUN = mean)
stat_DF[1, "CPM"] / stat_DF[2, "CPM"]
aggregate(CPM ~ group, data = expr_sub, FUN = function(x) { y <- sum(x > 0) / length(x) })

gp1 <- ggplot(expr_sub, aes(x = group, y = CPM, fill = group)) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "bar", show.legend = F) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) + 
  geom_path(data = data.frame(group = rep(c("Control", paste0("si", TF_id)), each = 2), CPM = 66 * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F), group = 1) + 
  annotate(geom = "text", x = 1.5, y = 66 * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  xlab(NULL) + ylab("CPM") + ggtitle(TF_id) + 
  scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(0, 66)) + 
  scale_fill_manual(values = c("grey", "#32CD32")) + 
  theme(aspect.ratio = 1, plot.title = element_text(face = "plain")) + theme(aspect.ratio = 2)

TF_id <- "NR2F1"
cell_sub <- subset(cellMetaData, group %in% c("Control", paste0("si", TF_id)), "cell", drop = T)
expr_sub <- as.data.frame(t(exprMatr_cellFiltered_CPM[TF_id, cell_sub, drop = F]))
expr_sub$group <- cellMetaData$group[match(rownames(expr_sub), cellMetaData$cell)]
colnames(expr_sub)[1] <- "CPM"
stat_DF <- aggregate(CPM ~ group, data = expr_sub, FUN = mean)
stat_DF[1, "CPM"] / stat_DF[2, "CPM"]
aggregate(CPM ~ group, data = expr_sub, FUN = function(x) { y <- sum(x > 0) / length(x) })

gp2 <- ggplot(expr_sub, aes(x = group, y = CPM, fill = group)) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "bar", show.legend = F) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) + 
  geom_path(data = data.frame(group = rep(c("Control", paste0("si", TF_id)), each = 2), CPM = 170 * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F), group = 1) + 
  annotate(geom = "text", x = 1.5, y = 170 * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  xlab(NULL) + ylab("CPM") + ggtitle(TF_id) + 
  scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(0, 170)) + 
  scale_fill_manual(values = c("grey", "#00B5EE")) + 
  theme(aspect.ratio = 1, plot.title = element_text(face = "plain")) + theme(aspect.ratio = 2)

pdf(file = paste0(OUT, "/TF_expr_inKD.pdf"), width = 5, height = 4, useDingbats = F)

plot_grid(gp1, gp2, align = "v")

dev.off()
