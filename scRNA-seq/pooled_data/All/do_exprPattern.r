# expression pattern
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

suppressMessages(library("arrow"))
library("ggplot2")
suppressMessages(library("cowplot"))
source("../../scripts/exprPattern_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/exprPattern/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

# load gene expression matrix (CPM)
expr_data_normed <- read_feather(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_cellFiltered_CPM.feather"))
expr_data_normed <- as.data.frame(expr_data_normed)
expr_data_normed_gene <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_cellFiltered_CPM.gene"), header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_data_normed) <- expr_data_normed_gene$V1

# load gene expression (avgCPM)
expr_data_avg <- read_feather(file = "03-expression/merged/filtering/UMIcount_cellFiltered_avgCPM_byCt.feather")
expr_data_avg <- as.data.frame(expr_data_avg)
expr_data_avg_gene <- read.table(file = "03-expression/merged/filtering/UMIcount_cellFiltered_avgCPM_byCt.gene", header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_data_avg) <- expr_data_avg_gene$V1

# load cell metatable
cellMetaData <- read.table(file = "cell_metatable_filtered_plus.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
dim(cellMetaData)
all(colnames(expr_data_normed) == rownames(cellMetaData))
cellMetaData$ts_ident <- paste(cellMetaData$tissue, cellMetaData$ident, sep = ".")

# load gene type
gene_type <- read.table("/rd/user/tianf/06-Human_cell_atlas/Genomes/human/gene_type_class.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 2)
colnames(gene_type) <- c("ensembl_id", "type", "class")
gene_type$class <- factor(gene_type$class, levels = unique(gene_type$class)[c(4,3,2,1,5)])

# load conservation (PhastCons and PhyloP)
gene_avgPc <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/UCSC/exon_avgPc.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_avgPc) <- c("ensembl_id", "gene", "avgPc")
gene_avgPp <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/UCSC/exon_avgPp.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_avgPp) <- c("ensembl_id", "gene", "avgPp")

# analysis
avgCPM_per_gene <- rowMeans(expr_data_avg)
nCellType_per_gene <- rowSums(expr_data_avg > 0)
pattern_stat <- data.frame(avgCPM = avgCPM_per_gene, nCellType = nCellType_per_gene, stringsAsFactors = F)
pattern_stat$class <- gene_type$class[match(rownames(pattern_stat), rownames(gene_type))]
pattern_stat$avgPc <- gene_avgPc$avgPc[match(rownames(pattern_stat), gene_avgPc$gene)]
pattern_stat$avgPp <- gene_avgPp$avgPp[match(rownames(pattern_stat), gene_avgPp$gene)]
pattern_stat$gene <- rownames(pattern_stat)
pattern_stat <- subset(pattern_stat, class %in% c("protein_coding", "lncRNA"))
pattern_stat$class <- factor(pattern_stat$class)
###
levels(pattern_stat$class)
levels(pattern_stat$class) <- Hmisc::capitalize(gsub("_", " ", levels(pattern_stat$class)))
levels(pattern_stat$class)
###

# sig text
wilcox.test(subset(pattern_stat, class == "Protein coding", "avgCPM", drop = T), 
            subset(pattern_stat, class == "LncRNA", "avgCPM", drop = T), alternative = "greater")$p.value

pdf(file = paste0(OUT, "/avgCPM_stat.pdf"), width = 5, height = 4, useDingbats = F)

ggplot(pattern_stat, aes(x = class, y = log10(avgCPM + 1))) + geom_violin(aes(fill = class, color = class), show.legend = F) + 
  stat_summary(fun.data = "mean_cl_boot", colour = "white") + 
  geom_line(data = data.frame(x = c(1,1,2,2), y = c(5,5.2,5.2,5)), aes(x = x, y = y)) + 
  annotate("text", x = 1.5, y = 5.4, label = "***", color = "red", size = 6) + 
  scale_fill_manual(values = c("hotpink", "dodgerblue")) + 
  scale_color_manual(values = c("hotpink", "dodgerblue")) + 
  xlab(NULL) + ylab(parse(text = "Log[10]~(CPM + 1)")) + 
  theme(aspect.ratio = 1)

ggplot(pattern_stat, aes(x = class, y = log10(avgCPM + 1))) + geom_violin(aes(fill = class, color = class), show.legend = F) + 
  stat_summary(geom = "point", fun = "median", colour = "white", size = 3) + 
  geom_line(data = data.frame(x = c(1, 1, 2, 2), y = c(5, 5.2, 5.2, 5)), aes(x = x, y = y)) + 
  annotate("text", x = 1.5, y = 5.4, label = "***", color = "red", size = 14) + 
  coord_cartesian(ylim = c(0, 5.6)) + 
  scale_fill_manual(values = c("hotpink", "dodgerblue")) + 
  scale_color_manual(values = c("hotpink", "dodgerblue")) + 
  xlab(NULL) + ylab(NULL) + 
  theme(aspect.ratio = 1, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 20))

# sort by avgCPM
pattern_stat_sorted <- pattern_stat[order(pattern_stat$avgCPM, decreasing = T), ]
pattern_stat_sorted$rank <- 1:nrow(pattern_stat_sorted)

pattern_stat_fmt <- pattern_stat_sorted
lncRNA_case <- head(subset(pattern_stat_fmt, class == "LncRNA", "gene", drop = T), 5)
pattern_stat_fmt_sub1 <- subset(pattern_stat_fmt, ! gene %in% lncRNA_case)
pattern_stat_fmt_sub1 <- pattern_stat_fmt_sub1[seq(1, nrow(pattern_stat_fmt_sub1), by = 8), ]
pattern_stat_fmt_sub2 <- subset(pattern_stat_fmt, gene %in% lncRNA_case)
pattern_stat_fmt_sub2$gene <- gsub("MALAT1--1", "MALAT1", pattern_stat_fmt_sub2$gene)
pattern_stat_fmt_sim <- rbind(pattern_stat_fmt_sub1, pattern_stat_fmt_sub2)

ggplot(pattern_stat_fmt_sim, aes(x = rank, y = log10(avgCPM + 1))) + geom_point(aes(color = class), alpha = 0.3) + 
  ggrepel::geom_text_repel(data = pattern_stat_fmt_sub2, aes(label = gene), force = 0.5, nudge_x = 2000, direction = "y", hjust = 0, point.padding = 0.5, min.segment.length = 0.1, color = "dodgerblue") + 
  scale_color_manual(values = c("hotpink", "dodgerblue")) + 
  xlab("Rank") + ylab(parse(text = "Log[10]~(CPM + 1)")) + labs(color = NULL) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  theme(aspect.ratio = 1, legend.position = c(0.495, 0.175), legend.justification = c(0, 0))

dev.off()

## add ligand/receptor
ligand <- read.table(file = "../../Data/CellPhoneDB/human_ligand.txt", header = F, stringsAsFactors = F)$V1
receptor <- read.table(file = "../../Data/CellPhoneDB/human_receptor.txt", header = F, stringsAsFactors = F)$V1

#pattern_stat_withLR <- pattern_stat
pattern_stat_ligand <- pattern_stat[rownames(pattern_stat) %in% ligand, ]
pattern_stat_ligand$class <- "Ligand"
pattern_stat_receptor <- pattern_stat[rownames(pattern_stat) %in% receptor, ]
pattern_stat_receptor$class <- "Receptor"
pattern_stat_withLR <- rbind(pattern_stat, rbind(pattern_stat_ligand, pattern_stat_receptor))
levels(pattern_stat_withLR$class)
levels(pattern_stat_withLR$class)[1] <- c("PCG")

pdf(file = paste0(OUT, "/avgCPM_stat_withLR.pdf"), width = 5, height = 4, useDingbats = F)

ggplot(pattern_stat_withLR, aes(x = class, y = log10(avgCPM + 1))) + geom_violin(aes(fill = class, color = class), show.legend = F) + 
  stat_summary(fun.data = "mean_cl_boot", colour = "white") + 
  #geom_line(data = data.frame(x = c(1,1,2,2), y = c(5,5.2,5.2,5)), aes(x = x, y = y)) + 
  #annotate("text", x = 1.5, y = 5.4, label = "***", color = "red", size = 6) + 
  #scale_fill_manual(values = c("hotpink", "dodgerblue")) + 
  #scale_color_manual(values = c("hotpink", "dodgerblue")) + 
  xlab(NULL) + ylab(parse(text = "Log[10]~(CPM + 1)")) + 
  theme(aspect.ratio = 1)

dev.off()
##

pdf(file = paste0(OUT, "/nCellType_stat.pdf"), width = 5, height = 4)

# ggplot(pattern_stat, aes(x = class, y = nCellType)) + geom_violin(aes(fill = class, color = class), show.legend = F) + 
#   stat_summary(fun.data = "mean_cl_boot", colour = "white") + 
#   scale_fill_manual(values = c("hotpink", "dodgerblue")) + 
#   scale_color_manual(values = c("hotpink", "dodgerblue")) + 
#   xlab(NULL) + ylab(parse(text = "Log[10]~(CPM + 1)")) + 
#   theme(aspect.ratio = 1)

ggplot(pattern_stat, aes(fill = class)) + 
  geom_histogram(data = subset(pattern_stat, class == "Protein coding"), aes(x = nCellType, y = ..density..), bins = 30, color = "white") + 
  geom_histogram(data = subset(pattern_stat, class == "LncRNA"), aes(x = nCellType, y = -..density..), bins = 30, color = "white") + 
  scale_fill_manual(values = c("hotpink", "dodgerblue"), drop = F) + 
  scale_x_continuous(expand = c(0.04, 0)) + 
  xlab("Number of cell type") + ylab("Density") + labs(fill = NULL) + 
  theme(aspect.ratio = 1, legend.position = c(0, 1.04), legend.justification = c(0, 1))

dev.off()

# gene_nCellType_ge150 <- subset(pattern_stat, nCellType >= 150, "gene", drop = T)
# write.table(x = gene_nCellType_ge150, file = file.path(OUT, "gene_nCellType_ge150.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
# gene_nCellType_ge200 <- subset(pattern_stat, nCellType >= 200, "gene", drop = T)
# write.table(x = gene_nCellType_ge200, file = file.path(OUT, "gene_nCellType_ge200.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
# gene_nCellType_ge228 <- subset(pattern_stat, nCellType >= 228, "gene", drop = T)
# write.table(x = gene_nCellType_ge228, file = file.path(OUT, "gene_nCellType_ge228.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
# gene_nCellType_le1 <- subset(pattern_stat, nCellType <= 1, "gene", drop = T)
# write.table(x = gene_nCellType_le1, file = file.path(OUT, "gene_nCellType_le1.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

pattern_stat_withCut <- pattern_stat
pattern_stat_withCut$nCellType_withCut <- cut(pattern_stat_withCut$nCellType, breaks = seq(0, 230, length.out = 6), include.lowest = T)
pattern_stat_withCut_rmNA <- subset(pattern_stat_withCut, ! is.na(avgPc))

# significance test
pvalue_LS <- apply(data.table::CJ(cut1 = levels(pattern_stat_withCut$nCellType_withCut), cut2 = levels(pattern_stat_withCut$nCellType_withCut), sorted = F), 1, function(i) {
  pvalue_pro_Pc <- wilcox.test(subset(pattern_stat_withCut, class == "Protein coding" & nCellType_withCut == i[1], "avgPc", drop = T), 
                               subset(pattern_stat_withCut, class == "Protein coding" & nCellType_withCut == i[2], "avgPc", drop = T), 
                               alternative = "less")$p.value
  pvalue_pro_Pp <- wilcox.test(subset(pattern_stat_withCut, class == "Protein coding" & nCellType_withCut == i[1], "avgPp", drop = T), 
                               subset(pattern_stat_withCut, class == "Protein coding" & nCellType_withCut == i[2], "avgPp", drop = T), 
                               alternative = "less")$p.value
  pvalue_lnc_Pc <- wilcox.test(subset(pattern_stat_withCut, class == "LncRNA" & nCellType_withCut == i[1], "avgPc", drop = T), 
                               subset(pattern_stat_withCut, class == "LncRNA" & nCellType_withCut == i[2], "avgPc", drop = T), 
                               alternative = "less")$p.value
  pvalue_lnc_Pp <- wilcox.test(subset(pattern_stat_withCut, class == "LncRNA" & nCellType_withCut == i[1], "avgPp", drop = T), 
                               subset(pattern_stat_withCut, class == "LncRNA" & nCellType_withCut == i[2], "avgPp", drop = T), 
                               alternative = "less")$p.value
  y <- data.frame(cut1 = i[1], cut2 = i[2], pvalue_pro_Pc, pvalue_pro_Pp, pvalue_lnc_Pc, pvalue_lnc_Pp, stringsAsFactors = F, row.names = NULL)
})
pvalue_DF <- do.call("rbind", pvalue_LS)
pvalue_DF$cut1 <- factor(pvalue_DF$cut1, levels = unique(pvalue_DF$cut1))
pvalue_DF$cut2 <- factor(pvalue_DF$cut2, levels = unique(pvalue_DF$cut2))
pheatmap::pheatmap((reshape2::acast(pvalue_DF, cut1 ~ cut2, value.var = "pvalue_pro_Pc") < 0.05) + 1, cluster_rows = F, cluster_cols = F, main = "Protein coding: Pc")
pheatmap::pheatmap((reshape2::acast(pvalue_DF, cut1 ~ cut2, value.var = "pvalue_pro_Pp") < 0.05) + 1, cluster_rows = F, cluster_cols = F, main = "Protein coding: Pp")
pheatmap::pheatmap((reshape2::acast(pvalue_DF, cut1 ~ cut2, value.var = "pvalue_lnc_Pc") < 0.05) + 1, cluster_rows = F, cluster_cols = F, main = "LncRNA: Pc")
pheatmap::pheatmap((reshape2::acast(pvalue_DF, cut1 ~ cut2, value.var = "pvalue_lnc_Pp") < 0.05) + 1, cluster_rows = F, cluster_cols = F, main = "LncRNA: Pp")
#

pdf(file = paste0(OUT, "/nCellType_Cons.pdf"), width = 5, height = 4)

x_labels <- gsub("-.*", "            ", gsub("\\[|\\]|\\(|\\)", "", gsub("," , "-", levels(pattern_stat_withCut_rmNA$nCellType_withCut))))
x_labels[5] <- gsub("      $", "228", x_labels[5])

ggplot(pattern_stat_withCut_rmNA, aes(x = nCellType_withCut, y = avgPc, fill = class, color = class)) + geom_split_violin(width = 1.4, show.legend = F) + 
  stat_summary(data = subset(pattern_stat_withCut_rmNA, class == "Protein coding"), geom = "point", 
               fun = "median", colour = "white", size = 3, position = position_nudge(x = -0.11), show.legend = F) + 
  stat_summary(data = subset(pattern_stat_withCut_rmNA, class == "LncRNA"), geom = "point", 
               fun = "median", colour = "white", size = 3, position = position_nudge(x = 0.082), show.legend = F) + 
  geom_line(data = data.frame(x = c(1, 1, 5, 5), y = c(1.1, 1.175, 1.175, 1.1)), aes(x = x, y = y), inherit.aes = F) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  #geom_vline(xintercept = seq(0.5, 5.5), linetype = "dashed") + 
  annotate("text", x = 3, y = 1.25, label = "***", color = "red", size = 6) + 
  scale_fill_manual(values = c("hotpink", "dodgerblue")) + 
  scale_color_manual(values = c("hotpink", "dodgerblue")) + 
  theme(aspect.ratio = 1, axis.ticks.x = element_line(color = NA)) + scale_x_discrete(labels = x_labels) + 
  labs(fill = NULL) + xlab("Number of cell type") + ylab("Average PhastCons in exon")

ggplot(pattern_stat_withCut_rmNA, aes(x = nCellType_withCut, y = avgPp, fill = class, color = class)) + geom_split_violin(width = 1.4, show.legend = F) + 
  stat_summary(data = subset(pattern_stat_withCut_rmNA, class == "Protein coding"), geom = "point", 
               fun = "median", colour = "white", size = 3, position = position_nudge(x = -0.11), show.legend = F) + 
  stat_summary(data = subset(pattern_stat_withCut_rmNA, class == "LncRNA"), geom = "point", 
               fun = "median", colour = "white", size = 3, position = position_nudge(x = 0.082), show.legend = F) + 
  geom_line(data = data.frame(x = c(1, 1, 5, 5), y = c(8, 8.5, 8.5, 8)), aes(x = x, y = y), inherit.aes = F) + 
  #scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  #geom_vline(xintercept = seq(0.5, 5.5), linetype = "dashed") + 
  annotate("text", x = 3, y = 9, label = "***", color = "red", size = 6) + 
  scale_fill_manual(values = c("hotpink", "dodgerblue")) + 
  scale_color_manual(values = c("hotpink", "dodgerblue")) + 
  theme(aspect.ratio = 1, axis.ticks.x = element_line(color = NA)) + scale_x_discrete(labels = x_labels) + 
  labs(fill = NULL) + xlab("Number of cell type") + ylab("Average PhyloP in exon")

dev.off()
