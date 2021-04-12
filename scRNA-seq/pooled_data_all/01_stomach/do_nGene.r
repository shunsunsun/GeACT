# detect gene number comparison
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/01_stomach/")

library("ggplot2")
library("cowplot")

samplingPos <- "."
OUT <- paste0("03-expression/merged/nGene/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

cell_metatable <- read.table("cell_metatable_filtered_plus.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cell_metatable)

cell_num_DF <- merge(as.data.frame(table(subset(cell_metatable, stage == "14w", "group", drop = T))), 
                     as.data.frame(table(subset(cell_metatable, stage == "20w", "group", drop = T))), by = "Var1")
colnames(cell_num_DF)[1] <- "group"
cell_metatable_sub <- subset(cell_metatable, group %in% subset(cell_num_DF, Freq.x >= 50 & Freq.y >= 50, "group", drop = T))
sapply(unique(cell_metatable_sub$group), function(i) {
  ka <- subset(cell_metatable, stage == "14w" & group == i, "nGene", drop = T)
  kb <- subset(cell_metatable, stage == "20w" & group == i, "nGene", drop = T)
  pv <- wilcox.test(x = ka, y = kb, alternative = "greater")$p.value
  return(pv)
})

gp <- ggplot(cell_metatable_sub, aes(x = stage, y = nGene, fill = stage)) + geom_violin() + facet_grid(. ~ group, switch = "x") + 
  annotate(geom = "segment", x = c(1,1,2), xend = c(1,2,2), y = c(7800,8000,8000), yend = c(8000,8000,7800)) + 
  annotate(geom = "text", x = 1.5, y = 8200, label = "***", color = "red", size = 5.5) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) + ylab("Detected gene number")
ggsave(filename = paste0(OUT, "/nGene.pdf"), plot = gp, width = 6, height = 5)
