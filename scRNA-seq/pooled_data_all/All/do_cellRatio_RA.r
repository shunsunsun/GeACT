
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("ggplot2")
library("cowplot")

cellMeta_RNA <- read.table(file = "cell_metatable_filtered_aligned.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_RNA)
cellMeta_ATAC <- read.table(file = "cell_metatable_ATAC.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_ATAC)

# load color
load("../../pooled_data/All/03-expression/merged/cellCluster/cg_color.RData")
#

cellRatio_RNA_LS <- lapply(unique(cellMeta_RNA$tissue), function(x) {
  tmp <- subset(cellMeta_RNA, tissue == x)
  y <- as.data.frame(table(tmp$group), stringsAsFactors = F)
  y$ratio <- y$Freq / sum(y$Freq)
  y$tissue <- x
  colnames(y)[1:2] <- c("group", "number")
  return(y)
})
cellRatio_RNA <- do.call("rbind", cellRatio_RNA_LS)

cellRatio_ATAC_LS <- lapply(unique(cellMeta_ATAC$tissue), function(x) {
  tmp <- subset(cellMeta_ATAC, tissue == x)
  y <- as.data.frame(table(tmp$group), stringsAsFactors = F)
  y$ratio <- y$Freq / sum(y$Freq)
  y$tissue <- x
  colnames(y)[1:2] <- c("group", "number")
  return(y)
})
cellRatio_ATAC <- do.call("rbind", cellRatio_ATAC_LS)
#
cellRatio_ATAC[cellRatio_ATAC$group == "TRUE", "group"] <- "T"
#

cellRatio_RA <- merge(cellRatio_RNA, cellRatio_ATAC, by = c("tissue", "group"), all = T, sort = F)
#
cellRatio_RA <- subset(cellRatio_RA, group != "Unknown")
cellRatio_RA <- subset(cellRatio_RA, tissue != "heart")
cellRatio_RA[is.na(cellRatio_RA$number.y), c("number.y", "ratio.y")] <- 0
#

cor_value <- cor(cellRatio_RA$ratio.x, cellRatio_RA$ratio.y)
cat("R^2:", cor_value^2, "\n")

cellRatio_RA_sub <- subset(cellRatio_RA, group %in% names(cg_color))
cellRatio_RA_sub$group <- factor(cellRatio_RA_sub$group, levels = names(cg_color))
cor_value <- cor(cellRatio_RA_sub$ratio.x, cellRatio_RA_sub$ratio.y)
cat("R^2:", cor_value^2, "\n")

pdf("cellRatio_RA.pdf", width = 6, height = 5)

ggplot(cellRatio_RA_sub, aes(x = ratio.x, y = ratio.y)) + geom_point(aes(color = group)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  annotate(geom = "text", x = 0, y = 1, label = paste0("R^2:~", round(cor_value^2, 2)), parse = T, hjust = 0) + 
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  xlab("RNA cell type proportion") + ylab("ATAC cell type proportion") + labs(color = NULL) + 
  scale_color_manual(values = cg_color) + 
  theme(aspect.ratio = 1)

dev.off()
