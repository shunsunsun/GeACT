# risk gene
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

suppressMessages(library("arrow"))
library("ggplot2")
suppressMessages(library("cowplot"))
suppressMessages(library("ComplexHeatmap"))

samplingPos <- "."
OUT <- paste0("03-expression/merged/riskGene/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/riskGene.RData"))

# 1. preprocess ----
# expression matrix
expr_data <- read_feather(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered.feather"))
dim(expr_data)
expr_data <- as.data.frame(expr_data)
expr_gene <- read.table(file = paste0("03-expression/merged/filtering/UMIcount_filtered.gene"), header = F, sep = "\t", stringsAsFactors = F, comment.char = "")
rownames(expr_data) <- expr_gene$V1

# cell meta
cellMetaData <- read.table(file = "cell_metatable_filtered_plus.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
dim(cellMetaData)
all(colnames(expr_CPM) == rownames(cellMetaData))
cellMetaData$ident.ori <- cellMetaData$ident
cellMetaData$ident <- paste(cellMetaData$tissue, cellMetaData$ident, sep = ".")

# cell type meta
ctMetaData <- read.table(file = "cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
dim(ctMetaData)
ctMetaData$ident.ori <- ctMetaData$ident
ctMetaData$ident <- paste(ctMetaData$tissue, ctMetaData$ident, sep = ".")
ctMetaData <- merge(ctMetaData, unique(cellMetaData[, c("ident", "group")]), by = "ident", sort = F)

# CPM
expr_CPM <- sweep(expr_data, MARGIN = 2, STATS = colSums(expr_data), FUN = "/") * 1e6
# avgCPM
expr_avgCPM <- sapply(split(rownames(cellMetaData), cellMetaData$ident), function(x) { y <- rowMeans(expr_CPM[, x, drop = F]); return(y) })
expr_aloCPM <- log10(expr_avgCPM + 1)

rm(expr_data, expr_gene); gc()

# read gene type info
gene_type <- read.table("/rd/user/tianf/06-Human_cell_atlas/Genomes/human/gene_type_class.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 2)
colnames(gene_type) <- c("ensembl_id", "type", "class")
gene_type$class <- factor(gene_type$class, levels = unique(gene_type$class)[c(4,3,2,1,5)])

# load conservation (PhastCons and PhyloP)
gene_avgPc <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/UCSC/exon_avgPc.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_avgPc) <- c("ensembl_id", "gene", "avgPc")
gene_avgPp <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/UCSC/exon_avgPp.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_avgPp) <- c("ensembl_id", "gene", "avgPp")

# risk genes
riskGene <- read.table(file = "../../Data/DisGeNET/riskGene_list.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(riskGene)
riskGene$diseaseName_new <- gsub(" of.*", "", riskGene$diseaseName)
riskGene[riskGene$tissue == "kidney" & riskGene$diseaseName_new == "Congenital anomaly", "diseaseName_new"] <- "Anomaly"
riskGene[riskGene$tissue == "kidney" & riskGene$diseaseName_new == "Unilateral agenesis", "diseaseName_new"] <- "Agenesis"
riskGene <- subset(riskGene, tissue != "liver")
riskGene <- subset(riskGene, ! geneSymbol %in% c("SUZ12", "PPP1R21"))

# 2. plot
pdf(paste0(OUT, "/riskGene_map.pdf"), width = 8, height = 6)

ht_opt(legend_title_gp = gpar(fontsize = 12), legend_labels_gp = gpar(fontsize = 12))
for(ts in unique(riskGene$tissue)) {
  cat(ts, "\n")
  riskGene_sub <- subset(riskGene, tissue == ts & geneSymbol %in% rownames(expr_aloCPM))
  cellType_sub <- subset(ctMetaData, tissue == ts)
  expr_sub <- expr_aloCPM[riskGene_sub$geneSymbol, gsub("\\..*", "", colnames(expr_aloCPM)) == ts, drop = F]
  expr_sub <- expr_sub[, match(cellType_sub$ident, colnames(expr_sub)), drop = F]
  colnames(expr_sub) <- gsub(".*\\.", "", colnames(expr_sub))

  ht <- Heatmap(matrix = t(expr_sub), name = "Expression", col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                border = "grey50", 
                cluster_rows = F, cluster_columns = T, show_parent_dend_line = F, 
                row_split = factor(cellType_sub$group, levels = unique(cellType_sub$group)), column_split = riskGene_sub$diseaseName_new, 
                row_title_gp = gpar(size = 12), row_title_rot = 0, 
                heatmap_legend_param = list(direction = "vertical", legend_height = unit(3, "cm")))
  print(ht)
}
ht_opt(RESET = T)

dev.off()

pdf(paste0(OUT, "/riskGene_map_diaphragm.pdf"), width = 6.5, height = 4)

ht_opt(legend_title_gp = gpar(fontsize = 12), legend_labels_gp = gpar(fontsize = 12))
for(ts in "diaphragm") {
  cat(ts, "\n")
  riskGene_sub <- subset(riskGene, tissue == ts & geneSymbol %in% rownames(expr_aloCPM))
  cellType_sub <- subset(ctMetaData, tissue == ts)
  expr_sub <- expr_aloCPM[riskGene_sub$geneSymbol, gsub("\\..*", "", colnames(expr_aloCPM)) == ts, drop = F]
  expr_sub <- expr_sub[, match(cellType_sub$ident, colnames(expr_sub)), drop = F]
  ###
  expr_sub <- expr_sub[- grep("^MYOD1$", rownames(expr_sub)), ]
  ###
  colnames(expr_sub) <- gsub(".*\\.", "", colnames(expr_sub))
  
  ht <- Heatmap(matrix = t(expr_sub), name = "Expression", col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                border = "grey50", 
                cluster_rows = F, cluster_columns = T, show_parent_dend_line = F, column_dend_height = unit(0.5, "cm"), 
                row_split = factor(cellType_sub$group, levels = unique(cellType_sub$group)), column_split = NULL, 
                row_title = NULL, row_title_gp = gpar(size = 12), row_title_rot = 0, show_heatmap_legend = F)
  lgd <- Legend(col_fun = circlize::colorRamp2(seq(min(expr_sub), max(expr_sub), length.out = 10), colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)), 
                title = "Expression", title_gp = gpar(fontsize = 12), title_gap = unit(3.5, "mm"), labels_gp = gpar(fontsize = 12), 
                at = seq(0, 3, by = 1), labels = seq(0, 3, by = 1), 
                direction = "vertical", legend_height = unit(3, "cm"))
  draw(ht, annotation_legend_list = list(lgd))
}
ht_opt(RESET = T)

dev.off()

expr_avgCPM_diaphragm <- expr_avgCPM[, grep("diaphragm", colnames(expr_avgCPM))]
gene_cor <- cor(t(expr_avgCPM_diaphragm))

# case gene
gene_case <- "GATA4"
gene_cor_case <- data.frame(case = gene_case, gene = colnames(gene_cor), cor_value = gene_cor[gene_case, ], stringsAsFactors = F)
gene_cor_case$class <- gene_type$class[match(gene_cor_case$gene, rownames(gene_cor_case))]
gene_cor_case <- gene_cor_case[order(gene_cor_case$cor_value, decreasing = T), ]
gene_cor_case_ftd <- subset(gene_cor_case, case != gene & (cor_value >= 0.957 | cor_value < -0.8))

gene_cor_sub <- gene_cor[c(gene_case, gene_cor_case_ftd$gene), c(gene_case, gene_cor_case_ftd$gene)]

gene_avgPc_case <- gene_avgPc$avgPc[match(rownames(gene_cor_sub), gene_avgPc$gene)]
gene_type_case <- gene_type$class[match(rownames(gene_cor_sub), rownames(gene_type))]

pdf(paste0(OUT, "/expr_cor_GATA4.pdf"), width = 5.5, height = 4)

ht <- Heatmap(matrix = gene_cor_sub, name = "Correlation", 
              col = circlize::colorRamp2(seq(-1, 1, length.out = 10), colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)), 
              border = "grey50", show_column_names = F, 
              cluster_rows = T, cluster_columns = T, show_parent_dend_line = F, row_dend_width = unit(0.5, "cm"), column_dend_height = unit(0.5, "cm"), 
              row_split = factor(1 - (rownames(gene_cor_sub) == gene_case)), column_split = factor(1 - (colnames(gene_cor_sub) == gene_case)), 
              row_title = NULL, column_title = NULL, show_heatmap_legend = F, 
              bottom_annotation = columnAnnotation(PhastCons = anno_barplot(gene_avgPc_case, height = unit(1.725, "cm"), border = F, 
                                                                            gp = gpar(fill = "#97de97", col = NA), 
                                                                            axis_param = list(at = c(0, 0.4, 0.8), labels = c(0, 0.4, 0.8), gp = gpar(fontsize = 12))), 
                                                   annotation_name_gp = gpar(fontsize = 12)))

lgd <- Legend(col_fun = circlize::colorRamp2(seq(-1, 1, length.out = 10), colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)), 
              title = "Correlation", title_gp = gpar(fontsize = 12), title_gap = unit(3.5, "mm"), labels_gp = gpar(fontsize = 12), 
              at = seq(-1, 1, by = 0.5), labels = seq(-1, 1, by = 0.5), 
              direction = "vertical", legend_height = unit(3, "cm"))
draw(ht, annotation_legend_list = list(lgd))

dev.off()

save.image(file = paste0(OUT, "/riskGene.RData"))
