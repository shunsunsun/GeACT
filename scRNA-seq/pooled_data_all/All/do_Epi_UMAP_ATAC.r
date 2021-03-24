# epithelial cells UMAP (ATAC)
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("ggplot2")
library("cowplot")

load("../../pooled_data/All/03-expression/merged/cellCluster/ct_color.RData")
tissue_ordered <- read.table(file = "../../pooled_data/All/tissue_ordered.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

meta <- read.table(file = "04-open_chromatin/merged/epi_dimred.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "", row.names = 1)
meta$tissue <- Hmisc::capitalize(meta$tissue)
meta$tissue <- factor(meta$tissue, levels = tissue_ordered)

pdf(file = "Epi_UMAP_ATAC.pdf", width = 6, height = 6)

ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = tissue)) + geom_point(size = 0.5, alpha = 0.6) + 
  scale_color_manual(values = ct_color[tissue_ordered]) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
  labs(color = NULL) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
  theme(aspect.ratio = 1)

dev.off()
