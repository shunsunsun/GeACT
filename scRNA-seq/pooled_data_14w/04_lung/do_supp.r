
TSNEPlot(object = expr_assigned, pt.size = 2, do.label = T, no.legend = T)

FeaturePlot(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "HBG1", "PECAM1", "COL1A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
# smooth muscle
FeaturePlot(object = expr, features.plot = c("CNN1", "ACTA2", "TAGLN", "NOTCH3", "LGR5", "DES", "LGR6", "MYH11", "PDGFRA"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# fibroblast
FeaturePlot(object = expr, features.plot = c("COL1A1", "PDGFRA", "ELN", "ACTA2", "PLIN2", "APOE"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# other
FeaturePlot(object = expr, features.plot = c("CSPG4", "TRPC6", "PDGFRB", "MSLN"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

FeaturePlot(object = expr, features.plot = c("CSPG4", "TRPC6", "PDGFRB", "MSLN"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

x <- "NOTCH3"
do_queryGene(x)
FeaturePlot(object = expr, features.plot = x, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.hover = T)
VlnPlot(object = expr, features.plot = c(x), size.title.use = 14, point.size.use = 0.1)

xxx <- FindMarkers(object = expr, ident.1 = 5, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx)

yyy <- FindMarkers(object = expr, ident.1 = 0, ident.2 = NULL, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(yyy)

zzz <- FindMarkers(object = expr, ident.1 = 2, ident.2 = 3, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(zzz)

kkk <- FindMarkers(object = expr, ident.1 = 12, ident.2 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(kkk)
kkk <- subset(kkk, power>=0.4 & avg_logFC>=log(2))
cat(rownames(kkk), sep = "\n")

kkk <- FindMarkers(object = expr, ident.1 = 10, ident.2 = 7, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
kkk <- subset(kkk, power>=0.4 & avg_logFC>=log(1.5))
cat(rownames(kkk), sep = "\n")

sss <- FindMarkers(object = expr, ident.1 = 14, ident.2 = c(7, 2, 1), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss)

xxx <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = 21, right1 = 26, top1 = 7, bottom1 = 2
                       )
head(xxx)

yyy <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = -38, right1 = -28, top1 = -7, bottom1 = -12.5, 
                       left2 = -38, right2= -28, top2 = -12.5, bottom2 = -18)
head(yyy)

tmp <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_unfiltered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
dim(tmp)
tmp <- tmp[, colnames(expr_data)]
xxx <- tmp[grep("^MUC", rownames(tmp)), ]
