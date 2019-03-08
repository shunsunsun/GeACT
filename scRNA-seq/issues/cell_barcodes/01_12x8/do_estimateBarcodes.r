
setwd("~/lustre/06-Human_cell_atlas/issues/barcodes")

# outer barcodes
outer_dist <- read.table("outer_barcodes_dist.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(outer_dist)
colnames(outer_dist) <- c("cmp1", "cmp2", "distance")
outer_dist$cmp1 <- factor(outer_dist$cmp1, levels = unique(outer_dist$cmp1))
outer_dist$cmp2 <- factor(outer_dist$cmp2, levels = unique(outer_dist$cmp2))

library("reshape2")
outer_dist_MT <- acast(data = outer_dist, formula = cmp1 ~ cmp2, value.var = "distance")
diag(outer_dist_MT) <- NA

pdf("outer_barcodes_dist.pdf", width = 6, height = 6, useDingbats = F, onefile = F)
library("pheatmap")
pheatmap(mat = outer_dist_MT, cluster_rows = F, cluster_cols = F, color = c("red", colorRampPalette(c("purple", "blue"))(5)), breaks = 1:7, display_numbers = T, number_format = "%d", fontsize_number = 14, number_color = "black")
dev.off()

# inner barcodes
inner_dist <- read.table("inner_barcodes_dist.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(inner_dist)
colnames(inner_dist) <- c("cmp1", "cmp2", "distance")
inner_dist$cmp1 <- factor(inner_dist$cmp1, levels = unique(inner_dist$cmp1))
inner_dist$cmp2 <- factor(inner_dist$cmp2, levels = unique(inner_dist$cmp2))

library("reshape2")
inner_dist_MT <- acast(data = inner_dist, formula = cmp1 ~ cmp2, value.var = "distance")
diag(inner_dist_MT) <- NA

pdf("inner_barcodes_dist.pdf", width = 6, height = 6, useDingbats = F, onefile = F)
library("pheatmap")
pheatmap(mat = inner_dist_MT, cluster_rows = F, cluster_cols = F, color = c("red", colorRampPalette(c("purple", "blue"))(5)), breaks = 1:7, display_numbers = T, number_format = "%d", fontsize_number = 14, number_color = "black")
dev.off()
