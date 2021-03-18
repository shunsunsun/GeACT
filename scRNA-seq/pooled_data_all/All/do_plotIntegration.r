
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("ggplot2")
library("cowplot")
library(grid)

# cell meta
cellMetaData <- read.table(file = "04-open_chromatin/merged/integration_dimreduc_19-22w.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
colnames(cellMetaData)[colnames(cellMetaData) == "celltype"] <- "ident"
cellMetaData$tech <- factor(cellMetaData$tech, levels = c("RNA", "ATAC"))

# cell type meta
ctMetaData <- read.table(file = "cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")

tissue_used <- intersect(unique(ctMetaData$tissue), unique(cellMetaData$tissue))
gp_LS <- lapply(tissue_used, function(ts_id) {
  cellMetaData_sub <- subset(cellMetaData, tissue == ts_id)
  ctMetaData_sub <- subset(ctMetaData, tissue == ts_id)
  point_size <- ifelse(nrow(cellMetaData_sub) >= 500, 0.8, 1.8)
  if(nrow(cellMetaData_sub) >= 10000) { point_size <- 0.5 }
  ident_labels <- paste(1:length(ctMetaData_sub$ident), ctMetaData_sub$ident, sep = ": ")
  gp_RA <- ggplot(cellMetaData_sub, aes(x = tSNE_1, y = tSNE_2, color = tech)) + geom_point(alpha = 0.6, size = point_size, show.legend = T) + 
    #scale_color_manual(labels = ident_labels, values = ctMetaData_sub$color) + 
    theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
    ggtitle(Hmisc::capitalize(ts_id)) + 
    theme(plot.title = element_text(face = "plain", size = 20)) + 
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
    labs(color = NULL)
  gp_ts <- ggplot(cellMetaData_sub, aes(x = tSNE_1, y = tSNE_2, color = ident)) + geom_point(alpha = 0.6, size = point_size, show.legend = T) + 
    scale_color_manual(labels = ident_labels, values = ctMetaData_sub$color) + 
    theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
    ggtitle(Hmisc::capitalize(ts_id)) + 
    theme(plot.title = element_text(face = "plain", size = 20)) + 
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
    labs(color = NULL)
  y <- list(gp_RA, gp_ts)
  return(y)
})
names(gp_LS) <- tissue_used

#plot_grid(plotlist = gp_LS, ncol = 3, align = "hv")
gp_RA_all <- lapply(gp_LS, function(x) { y <- x[[1]] + theme(legend.position = "none", plot.title = element_text(face = "plain", size = 16))})
gp_ts_all <- lapply(gp_LS, function(x) { y <- x[[2]] + theme(legend.position = "none", plot.title = element_text(face = "plain", size = 16))})

pdf("integration_tSNE_byTissue.pdf", width = 21, height = 29.7, useDingbats = F)
plot_grid(plotlist = gp_RA_all, ncol = 3, align = "hv")
dev.off()

pdf("integration_tSNE_small_intestine.pdf", width = 8, height = 4, useDingbats = F)
plot_grid(plotlist = list(gp_RA_all[["small intestine"]], gp_ts_all[["small intestine"]]), ncol = 2, align = "hv")
dev.off()

pdf("integration_tSNE_4ts.pdf", width = 8, height = 8, useDingbats = F)

plot_grid(plotlist = gp_RA_all[c("small intestine", "pancreas", "kidney", "ovary")], ncol = 2, align = "hv")

# legend
lg <- cowplot::get_legend(gp_LS[[1]][[1]])
grid.newpage()
grid.draw(lg)

dev.off()
