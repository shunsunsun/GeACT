# plot combined tSNE or UMAP
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("ggplot2")
suppressMessages(library("cowplot"))

# cell meta
cellMetaData <- read.table(file = "cell_metatable_RNA_global.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)

# cell type meta
ctMetaData <- read.table(file = "cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")

gp_LS <- lapply(unique(ctMetaData$tissue), function(ts_id) {
  cellMetaData_sub <- subset(cellMetaData, tissue == ts_id)
  ctMetaData_sub <- subset(ctMetaData, tissue == ts_id)
  cellMetaData_sub$ident <- factor(cellMetaData_sub$ident, levels = ctMetaData_sub$ident)
  point_size <- ifelse(nrow(cellMetaData_sub) >= 500, 0.8, 1.8)
  if(nrow(cellMetaData_sub) >= 10000) { point_size <- 0.5 }
  ident_labels <- paste(1:length(ctMetaData_sub$ident), ctMetaData_sub$ident, sep = ": ")
  gp_ts <- ggplot(cellMetaData_sub, aes(x = tSNE_1_ali, y = tSNE_2_ali, color = ident)) + geom_point(alpha = 0.6, size = point_size) + 
    scale_color_manual(labels = ident_labels, values = ctMetaData_sub$color) + 
    theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
    ggtitle(Hmisc::capitalize(ts_id)) + 
    theme(plot.title = element_text(face = "plain", size = 20)) + 
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) + 
    labs(color = NULL) 
  return(gp_ts)
})

pdf("tSNE_byTissue.pdf", width = 21, height = 29.7, useDingbats = F)
plot_grid(plotlist = gp_LS, ncol = 3, align = "hv")
dev.off()
