
TSNEPlot(object = expr_assigned, pt.size = 2, do.label = T, no.legend = T)

mk_20w <- read.table(file = paste0(gsub("_14w", "", getwd()), "/03-expression/merged/cellCluster/Seurat_markerGenes.txt"), header = T, sep = "\t", stringsAsFactors = F)
dim(mk_20w)
mk_20w_top1 <- mk_20w %>% group_by(ident) %>% top_n(1, avg_logFC)

pdf(paste0(OUT, "/Seurat_mk_20w.pdf"), width = 6, height = 6, useDingbats = F)
for(i in mk_20w_top1$gene) {
  FeaturePlot(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
}
dev.off()

xxx <- read.table("03-expression/merged/filtering/filtering_cells.txt", header = T, sep = "\t", stringsAsFactors = F)
ggplot(subset(xxx, filter), aes(x = nUMI, y = nGene, color = X %in% "XCX_C-A6_01")) + geom_point(alpha = 0.6)

FeaturePlot(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "HBG1", "PECAM1", "COL1A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
# smooth muscle
FeaturePlot(object = expr, features.plot = c("CNN1", "ACTA2", "TAGLN", "NOTCH3", "LGR5", "DES", "LGR6", "MYH11"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# fibroblast
FeaturePlot(object = expr, features.plot = c("COL1A1", "PDGFRA", "ELN", "ACTA2", "PLIN2", "APOE"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# other
FeaturePlot(object = expr, features.plot = c("PLP1", "S100B", "HBG1","CSPG4", "TRPC6", "PDGFRB"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

FeaturePlot(object = expr, features.plot = c("PLP1", "S100B"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

FeaturePlot(object = expr, features.plot = c("HBG1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

x <- c("TOP2A")
do_queryGene(x)
#x %in% expr.markers$gene
#expr.markers[grep(x, expr.markers$gene), ]
FeaturePlot(object = expr, features.plot = x, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.hover = T)
VlnPlot(object = expr, features.plot = c(x), size.title.use = 14, point.size.use = 0.1)

# c small
xxx <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = 13, right1 = 16, bottom1 = 26, top1 = 29)
head(xxx)

xxx1 <- FindMarkers(object = expr, ident.1 = 11, ident.2 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx1)

xxx2 <- FindMarkers(object = expr, ident.1 = c(2,3,6), ident.2 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx2)

yyy <- FindMarkers(object = expr, ident.1 = 19, ident.2 = 11, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(yyy)

zzz <- FindMarkers(object = expr, ident.1 = 7, ident.2 = 15, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(zzz)

kkk <- FindMarkers(object = expr, ident.1 = 12, ident.2 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(kkk)
kkk <- subset(kkk, power>=0.4 & avg_logFC>=log(2))
cat(rownames(kkk), sep = "\n")

kkk <- FindMarkers(object = expr, ident.1 = 14, ident.2 = 10, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(kkk)
kkk <- subset(kkk, power>=0.4 & avg_logFC>=log(1.5))
cat(rownames(kkk), sep = "\n")

sss1 <- FindMarkers(object = expr, ident.1 = 0, ident.2 = 2, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss1)

sss2 <- FindMarkers(object = expr, ident.1 = 2, ident.2 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss2)

xxx <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = 5.9, right1 = 18, top1 = 12.5, bottom1 = 3, 
                       left2 = NULL, right2= NULL, top2 = NULL, bottom2 = NULL
                       )
head(xxx)

zzz <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = -17, right1 = -11, top1 = 33, bottom1 = 26, 
                       left2 = -10, right2 = -3.9, top2 = 32, bottom2 = 20.5)
head(zzz)

kkk <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left2 = -17, right2 = -11, top2 = 33, bottom2 = 29, 
                       left1 = -17, right1 = -11, top1 = 28, bottom1 = 26)
head(kkk)

tmp <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_unfiltered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
dim(tmp)
tmp <- tmp[, colnames(expr_data)]
xxx <- tmp[grep("^MUC", rownames(tmp)), ]
