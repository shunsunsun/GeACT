# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

#library("Seurat")
#library("dplyr")
library("Matrix")
library("pheatmap")
library("reshape2")
library("grid")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
#library("topGO")
library("RColorBrewer")
source("../../scripts/cluster_tools.r")
source("../../scripts/pheatmap_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/cellCluster/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

cellMetaData <- read.table("cell_metatable_ATAC.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMetaData)
cellMetaData$tissue <- Hmisc::capitalize(cellMetaData$tissue)
cellMetaData$group_new <- cellMetaData$group
#cellMetaData$group_new <- gsub("lium$", "lial", cellMetaData$group_new)
group_shown <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "B", 
                 "DC/Macrophage", "NKT", "T", "Glial", "FGC", 
                 "Granulosa", "Sertoli", "Erythrocyte", "Other")
cellMetaData$group_new[! cellMetaData$group_new %in% group_shown] <- "Other"
cellMetaData$group_new <- factor(cellMetaData$group_new, levels = group_shown)

load(paste0(OUT, "/ct_color.RData"))
ct_color <- ct_color[unique(cellMetaData$tissue)]

load(paste0(OUT, "/cg_color.RData"))
cg_color <- cg_color[names(cg_color) %in% cellMetaData$group_new]

pdf(paste0(OUT, "/Seurat_tSNE_ATAC.pdf"), width = 6, height = 6, useDingbats = F)

# color by tissue
ggplot(cellMetaData, aes(x = all_tSNE_1, y = all_tSNE_2, color = tissue)) + geom_point(size = 0.2) +
  theme(legend.position = "bottom", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, legend.box.margin = margin(t = -15), plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color)

dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_legend_1col_ATAC.pdf"), width = 9, height = 6, useDingbats = F)

# color by tissue
ggplot(cellMetaData, aes(x = all_tSNE_1, y = all_tSNE_2, color = tissue)) + geom_point(size = 0.2) +
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.8, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color)

dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_and_legend_ATAC.pdf"), width = 6, height = 5, useDingbats = F)

# color by tissue
gp <- ggplot(cellMetaData, aes(x = all_tSNE_1, y = all_tSNE_2, color = tissue)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.725, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color) + 
  ggtitle("Chromatin accessibility") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

# color by cell group
gp <- ggplot(cellMetaData, aes(x = all_tSNE_1, y = all_tSNE_2, color = group_new)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.8, "cm")) + 
  labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = cg_color) + 
  ggtitle("Chromatin accessibility") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

dev.off()

pdf(paste0(OUT, "/Seurat_UMAP_and_legend_ATAC.pdf"), width = 6, height = 5, useDingbats = F)

# color by tissue
gp <- ggplot(cellMetaData, aes(x = all_UMAP_1, y = all_UMAP_2, color = tissue)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.725, "cm")) + 
  labs(color = NULL) + xlab("UMAP1") + ylab("UMAP2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = ct_color) + 
  ggtitle("Chromatin accessibility") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

# color by cell group
gp <- ggplot(cellMetaData, aes(x = all_UMAP_1, y = all_UMAP_2, color = group_new)) + geom_point(size = 0.2) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.justification = c("center")) + 
  theme(aspect.ratio = 1, plot.margin = margin(t = 10, l = -20, r = 20)) + 
  theme(legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.8, "cm")) + 
  labs(color = NULL) + xlab("UMAP1") + ylab("UMAP2") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))) + 
  scale_color_manual(values = cg_color) + 
  ggtitle("Chromatin accessibility") + theme(plot.title = element_text(face = "plain"))
gp + theme(legend.position = "none")

# only legend
lg <- cowplot::get_legend(gp)
grid.newpage()
grid.draw(lg)

dev.off()

pdf(paste0(OUT, "/Seurat_tSNE_simple_ATAC.pdf"), width = 6, height = 6, useDingbats = F)

# color by tissue
ggplot(cellMetaData, aes(x = all_tSNE_1, y = all_tSNE_2, color = tissue)) + geom_point(size = 0.2, show.legend = F) +
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

png(paste0(OUT, "/Seurat_tSNE_simple_ATAC.png"), bg = "transparent", width = 3000, height = 3000, res = 600)

# color by tissue
ggplot(cellMetaData, aes(x = all_tSNE_1, y = all_tSNE_2, color = tissue)) + geom_point(size = 0.2, show.legend = F) +
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