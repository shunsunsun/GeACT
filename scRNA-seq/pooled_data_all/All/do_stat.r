
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("ggplot2")
suppressMessages(library("cowplot"))
library("ggridges")
library("ggsci")

load("../../pooled_data/All/03-expression/merged/cellCluster/ct_color.RData")

### RNA-seq ----
cellMeta <- read.table("cell_metatable.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta)
cellMeta$tissue <- gsub("_", " ", Hmisc::capitalize(cellMeta$tissue))
cellMeta$tissue <- factor(cellMeta$tissue, levels = sort(unique(cellMeta$tissue)))
cellMeta_ftd <- subset(cellMeta, QC)
ts_list1 <- levels(cellMeta_ftd$tissue)
length(ts_list1)
cat("QC ratio:", paste0(round(nrow(cellMeta_ftd)/nrow(cellMeta)*100, 2), "%"), paste0("(", nrow(cellMeta_ftd), "/", nrow(cellMeta), ")"), "\n")

# attr stat
cat("Clean reads (median):", median(cellMeta_ftd$cleanReads) / 1e6, "million", "\n")
cat("Genes (median):", median(cellMeta_ftd$nGene), "\n")
cat("UMI (median):", median(cellMeta_ftd$nUMI), "\n")

cellMeta_tb <- as.data.frame(table(cellMeta_ftd$tissue))
colnames(cellMeta_tb) <- c("tissue", "cellNum")

pdf("cellNum_stat.pdf", width = 4, height = 4)

# ggplot(cellMeta_tb, aes(x = tissue, y = log10(cellNum), fill = tissue)) + geom_bar(stat = "identity", show.legend = F) + coord_flip() + 
#   scale_x_discrete(limits = rev(levels(cellMeta_ftd$tissue))) + 
#   scale_y_continuous(limits = c(0, max(log10(cellMeta_tb$cellNum)) * 1.05), expand = c(0, 0)) + 
#   xlab("Organ") + ylab(expression(Log[10] ~ "(Cell number)")) + 
#   theme(plot.subtitle = element_text(hjust = 0.5)) #+ 
#   #ggtitle(label = "scRNA-seq", subtitle = paste("N =", format(nrow(cellMeta_ftd), big.mark = ",", scientific = F)))

ggplot(cellMeta_tb, aes(x = tissue, y = log10(cellNum), fill = tissue)) + geom_bar(stat = "identity", alpha = 0.6, show.legend = F) + coord_flip() + 
  scale_x_discrete(limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_y_continuous(limits = c(0, max(log10(cellMeta_tb$cellNum)) * 1.05), expand = c(0, 0), 
                     breaks = 0:4, labels = parse(text = paste0(10, "^", 0:4))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_tb$tissue)]) + 
  xlab("Organ") + ylab("Cell number") + 
  theme(plot.subtitle = element_text(hjust = 0.5))

dev.off()

# cellMeta_melted <- reshape2::melt(table(cellMeta[, c("tissue", "QC")]))
# ggplot(cellMeta_melted, aes(x = tissue, y = value, fill = tissue, alpha = QC)) + geom_bar(stat = "identity", show.legend = F) + 
#   scale_alpha_discrete(range = c(0.4,1)) + 
#   scale_y_continuous(limits = c(0, max(table(cellMeta$tissue)) * 1.05), expand = c(0, 0)) + xlab("Organ") + ylab("Cell number") + 
#   geom_text(aes(label = value, y = value), size = 3, position = position_stack(vjust = 0.5), show.legend = F) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("expr_stat.pdf", width = 4, height = 4)

# ggplot(cellMeta_ftd, aes(x = tissue, y = nGene, fill = tissue)) + geom_violin(show.legend = F) + coord_flip() + 
#   scale_x_discrete(limits = rev(levels(cellMeta_ftd$tissue))) + 
#   xlab("Organ") + ylab("Cell number") + 
#   theme(plot.margin = margin(r = 18), aspect.ratio = 0.5) 
# 
# ggplot(cellMeta_ftd, aes(x = nGene, fill = tissue)) + geom_density(show.legend = F) + facet_grid(tissue ~ .) + 
#   xlab("Organ") + ylab("Cell number")

ggplot(cellMeta_ftd, aes(x = log10(cleanReads), y = tissue, fill = tissue, color = tissue)) + geom_density_ridges(alpha = 0.6, show.legend = F) + 
  scale_x_continuous(breaks = 0:7, labels = parse(text = paste0(10, "^", 0:7))) + 
  scale_y_discrete(expand = c(0.01, 0), limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  scale_color_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  xlab("Clean reads per cell") + ylab("Organ")

ggplot(cellMeta_ftd, aes(x = nGene, y = tissue, fill = tissue, color = tissue)) + geom_density_ridges(alpha = 0.6, show.legend = F) + 
  scale_y_discrete(expand = c(0.01, 0), limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  scale_color_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  xlab("Genes per cell") + ylab("Organ") + 
  coord_cartesian(xlim = c(0, 1e4)) + theme(plot.margin = margin(7,10,7,7))

ggplot(cellMeta_ftd, aes(x = log10(nUMI), y = tissue, fill = tissue, color = tissue)) + geom_density_ridges(alpha = 0.6, show.legend = F) + 
  scale_x_continuous(breaks = 0:7, labels = parse(text = paste0(10, "^", 0:7))) + 
  scale_y_discrete(expand = c(0.01, 0), limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  scale_color_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  xlab("UMIs per cell") + ylab("Organ")

dev.off()

### ATAC-seq ----
cellMeta <- read.table("cell_metatable_ATAC.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta)
cellMeta <- cellMeta[, c("tissue", "cleanReads", "nFrags", "nPeak")]
colnames(cellMeta) <- c("tissue", "cleanReads", "frag_num", "peak_num")
cellMeta$tissue <- gsub("_", " ", Hmisc::capitalize(cellMeta$tissue))
cellMeta$tissue <- factor(cellMeta$tissue, levels = sort(unique(cellMeta$tissue)))
cellMeta_ftd <- cellMeta
ts_list2 <- levels(cellMeta_ftd$tissue)
length(ts_list2)
#cat("QC ratio:", paste0(round(nrow(cellMeta_ftd)/nrow(cellMeta)*100, 2), "%"), paste0("(", nrow(cellMeta_ftd), "/", nrow(cellMeta), ")"), "\n")

# attr stat
cat("Clean reads (median):", round(median(cellMeta_ftd$cleanReads)), "\n")
cat("Fragments (median):", median(cellMeta_ftd$frag_num), "\n")
cat("Peaks (median):", median(cellMeta_ftd$peak_num), "\n")

cellMeta_tb <- as.data.frame(table(cellMeta_ftd$tissue))
colnames(cellMeta_tb) <- c("tissue", "cellNum")
#cellMeta_tb$tissue <- factor(cellMeta_tb$tissue, levels = ts_list1)

pdf("cellNum_stat_ATAC.pdf", width = 4, height = 4)

# ggplot(cellMeta_tb, aes(x = tissue, y = log10(cellNum), fill = tissue)) + geom_bar(stat = "identity", show.legend = F) + coord_flip() + 
#   scale_fill_discrete(drop = F) + 
#   scale_x_discrete(limits = rev(levels(cellMeta_ftd$tissue))) + 
#   scale_y_continuous(limits = c(0, max(log10(cellMeta_tb$cellNum)) * 1.05), expand = c(0, 0)) + 
#   xlab("Organ") + ylab(expression(Log[10] ~ "(Cell number)")) + 
#   theme(plot.subtitle = element_text(hjust = 0.5)) #+ 
#   #ggtitle(label = "scRNA-seq", subtitle = paste("N =", format(nrow(cellMeta_ftd), big.mark = ",", scientific = F)))

ggplot(cellMeta_tb, aes(x = tissue, y = log10(cellNum), fill = tissue)) + geom_bar(stat = "identity", alpha = 0.6, show.legend = F) + coord_flip() + 
  scale_x_discrete(limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_y_continuous(limits = c(0, max(log10(cellMeta_tb$cellNum)) * 1.15), expand = c(0, 0), 
                     breaks = 0:4, labels = parse(text = paste0(10, "^", 0:4))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_tb$tissue)]) + 
  xlab("Organ") + ylab("Cell number") + 
  theme(plot.subtitle = element_text(hjust = 0.5))

dev.off()

pdf("expr_stat_ATAC.pdf", width = 4, height = 4)

ggplot(cellMeta_ftd, aes(x = log10(frag_num), y = tissue, fill = tissue, color = tissue)) + geom_density_ridges(alpha = 0.6, show.legend = F) + 
  scale_x_continuous(breaks = 0:7, labels = parse(text = paste0(10, "^", 0:7)), limits = c(3.8, 6)) + 
  scale_y_discrete(expand = c(0.01, 0), limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  scale_color_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  xlab("Fragments per cell") + ylab("Organ")

ggplot(cellMeta_ftd, aes(x = log10(peak_num), y = tissue, fill = tissue, color = tissue)) + geom_density_ridges(alpha = 0.6, show.legend = F) +
  scale_x_continuous(breaks = 0:7, labels = parse(text = paste0(10, "^", 0:7)), limits = c(3, 5)) + 
  scale_y_discrete(expand = c(0.01, 0), limits = rev(levels(cellMeta_ftd$tissue))) + 
  scale_fill_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  scale_color_manual(values = ct_color[levels(cellMeta_ftd$tissue)]) + 
  xlab("Peaks per cell") + ylab("Organ")

dev.off()
