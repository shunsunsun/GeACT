# plot combined tSNE or UMAP
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/01_stomach/")

library("ggplot2")
suppressMessages(library("cowplot"))

samplingPos <- "."
OUT <- paste0("03-expression/merged/dimReduction/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

# cell meta
cellMetaData <- read.table(file = "cell_metatable_filtered_plus.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
cellMetaData$samplingPos <- Hmisc::capitalize(cellMetaData$samplingPos)

# cell type meta
ctMetaData <- read.table(file = "cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
ctMetaData <- subset(ctMetaData, stage == "14w")[-4]
ctMetaData <- merge(ctMetaData, unique(cellMetaData[, c("ident", "group")]), by = "ident", sort = F)

# CellBlast result
ts_id <- Hmisc::capitalize(gsub(".*[0-9][0-9]_", "", getwd()))
cblast_res <- read.table(file = paste0("../../pooled_data_14w/cblast_result/", ts_id, ".csv"), header = T, sep = ",", stringsAsFactors = F, row.names = 1)
rownames(cblast_res) <- gsub(".", "-", rownames(cblast_res), fixed = T)

# merge info
cellMetaData <- merge(cellMetaData, cblast_res[, c("tSNE1", "tSNE2", "UMAP1", "UMAP2")], by = 0, sort = F)
colnames(cellMetaData)[c(1, 23:26)] <- c("cell", "tSNE_1M", "tSNE_2M", "UMAP_1M", "UMAP_2M")

cellMetaData$samplingPos <- factor(cellMetaData$samplingPos, levels = c("Fundus", "Body", "Antrum"))
cellMetaData$ident <- factor(cellMetaData$ident, levels = ctMetaData$ident)
cellMetaData$group <- factor(cellMetaData$group, levels = unique(ctMetaData$group))

ident_labels <- paste(1:length(levels(cellMetaData$ident)), levels(cellMetaData$ident), sep = ": ")

pdf(file = file.path(OUT, "tSNE_UMAP.pdf"), width = 9, height = 5.5, useDingbats = F)

# tSNE
ggplot(cellMetaData, aes(x = tSNE_1M, y = tSNE_2M, color = ident)) + geom_point(alpha = 0.6, size = 0.5) + 
  scale_color_manual(labels = ident_labels, values = ctMetaData$color) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + xlab("tSNE-1") + ylab("tSNE-2")

ggplot(cellMetaData, aes(x = tSNE_1M, y = tSNE_2M, color = stage)) + geom_point(alpha = 0.6, size = 0.5) + 
  #scale_color_manual(labels = ident_labels, values = ctMetaData$color) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + xlab("tSNE-1") + ylab("tSNE-2")

ggplot(cellMetaData, aes(x = tSNE_1M, y = tSNE_2M, color = samplingPos)) + geom_point(alpha = 0.6, size = 0.5) + 
  #scale_color_manual(labels = ident_labels, values = ctMetaData$color) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + xlab("tSNE-1") + ylab("tSNE-2")

# UMAP
ggplot(cellMetaData, aes(x = UMAP_1M, y = UMAP_2M, color = ident)) + geom_point(alpha = 0.6, size = 0.5) + 
  scale_color_manual(labels = ident_labels, values = ctMetaData$color) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + xlab("UMAP-1") + ylab("UMAP-2")

ggplot(cellMetaData, aes(x = UMAP_1M, y = UMAP_2M, color = stage)) + geom_point(alpha = 0.6, size = 0.5) + 
  #scale_color_manual(labels = ident_labels, values = ctMetaData$color) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + xlab("UMAP-1") + ylab("UMAP-2")

ggplot(cellMetaData, aes(x = UMAP_1M, y = UMAP_2M, color = samplingPos)) + geom_point(alpha = 0.6, size = 0.5) + 
  #scale_color_manual(labels = ident_labels, values = ctMetaData$color) + 
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), panel.border = element_blank()) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + xlab("UMAP-1") + ylab("UMAP-2")

dev.off()
