
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)

mk_20w <- read.table(file = paste0(gsub("_14w", "", getwd()), "/03-expression/merged/cellCluster/Seurat_markerGenes.txt"), header = T, sep = "\t", stringsAsFactors = F)
dim(mk_20w)
mk_20w_top1 <- mk_20w %>% group_by(ident) %>% top_n(1, avg_logFC)

pdf(paste0(OUT, "/Seurat_mk_20w.pdf"), width = 6, height = 6, useDingbats = F)
for(i in mk_20w_top1$gene) {
  FeaturePlot(object = expr, features.plot = i, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
}
dev.off()

FeaturePlot(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "HBG1", "PECAM1", "COL1A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
# smooth muscle
FeaturePlot(object = expr, features.plot = c("CNN1", "ACTA2", "TAGLN", "NOTCH3", "LGR5", "DES", "LGR6", "MYH11"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# fibroblast
FeaturePlot(object = expr, features.plot = c("COL1A1", "PDGFRA", "ELN", "ACTA2", "PLIN2", "APOE"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# other
FeaturePlot(object = expr, features.plot = c("CSPG4", "TRPC6", "PDGFRB"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

x <- "KIT"
do_queryGene(x)
FeaturePlot(object = expr, features.plot = x, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.hover = F)
VlnPlot(object = expr, features.plot = c(x), size.title.use = 14, point.size.use = 0.1)

# c3
xxx1 <- FindMarkers(object = expr, ident.1 = 3, ident.2 = NULL, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx1)
# c4
xxx2 <- FindMarkers(object = expr, ident.1 = 4, ident.2 = NULL, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx2)
xxy2 <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                        left1 = -8, right1 = 2, bottom1 = 20, top1 = 30, 
                        left2 = 2, right2 = 15, bottom2 = 23, top2 = 40)
head(xxy2)
xxz2 <- FindMarkers(object = expr, ident.1 = 15, ident.2 = 1, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxz2)
xxz2 <- FindMarkers(object = expr, ident.1 = 15, ident.2 = 4, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxz2)
xxz2 <- FindMarkers(object = expr, ident.1 = 15, ident.2 = c(1,4), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxz2)
FeaturePlot(object = expr, features.plot = rownames(xxz2)[1:4], cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.hover = F)
# c7
xxz3 <- FindMarkers(object = expr, ident.1 = 16, ident.2 = 7, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxz3)
# c8
xxz4 <- FindMarkers(object = expr, ident.1 = 17, ident.2 = 8, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxz4)
xxz4 <- FindMarkers(object = expr, ident.1 = 8, ident.2 = 17, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxz4)
# fibro
yyy <- FindMarkers(object = expr, ident.1 = c(1,9,15), ident.2 = c(0,2,5), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
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
                       left1 = 10, right1 = 20, top1 = -35, bottom1 = -40, 
                       left2 = NULL, right2= NULL, top2 = NULL, bottom2 = NULL
                       )
head(xxx)

zzz <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = -14, right1 = -12, bottom1 = 26, top1 = 29)
head(zzz)

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
