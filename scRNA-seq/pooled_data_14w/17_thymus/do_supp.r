
TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)

FeaturePlot(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "PECAM1", "COL1A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
# smooth muscle
FeaturePlot(object = expr, features.plot = intersect(c("CNN1", "ACTA2", "TAGLN", "NOTCH3", "LGR5", "DES", "LGR6", "MYH11", "PDGFRA"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# fibroblast
FeaturePlot(object = expr, features.plot = intersect(c("COL1A1", "PDGFRA", "ELN", "ACTA2", "PLIN2", "APOE"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# other
FeaturePlot(object = expr, features.plot = intersect(c("CSPG4", "TRPC6", "PDGFRB", "MSLN"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# RB
FeaturePlot(object = expr, features.plot = intersect(c("HBG1", "HBB"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

## https://www.abcam.com/primary-antibodies/t-cells-basic-immunophenotyping
# Th1
FeaturePlot(object = expr, features.plot = intersect(c("CXCR3", "IFNG", "IL2", "IL12A", "IL18", "STAT4", "STAT1"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Th2
FeaturePlot(object = expr, features.plot = intersect(c("CCR4", "IL2", "IL4"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Treg
FeaturePlot(object = expr, features.plot = intersect(c("CD4", "IL2RA", "IL7R", "CTLA4", "FOXP3", "STAT5A"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Killer T cell
FeaturePlot(object = expr, features.plot = intersect(c("CD8A", "EOMES", "TNF"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
##


# B
FeaturePlot(object = expr, features.plot = intersect(c("CD79A", "CD24", "MS4A1", "CD19"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Plasma
FeaturePlot(object = expr, features.plot = intersect(c("CD79A", "CD27", "SLAMF7"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# T naive
FeaturePlot(object = expr, features.plot = intersect(c("CD3E", "CD4", "CD8A", "GZMH", "GZMB"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# T meme/eff
FeaturePlot(object = expr, features.plot = intersect(c("COTL1", "LDHB"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# T 
FeaturePlot(object = expr, features.plot = intersect(c("CCR7", "LEF1"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# NK
FeaturePlot(object = expr, features.plot = intersect(c("KLRD1", "NKG7", "TYROBP", "GNLY"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# NK T
FeaturePlot(object = expr, features.plot = intersect(c("FCER1G", "TYROBP"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Neutrophil
FeaturePlot(object = expr, features.plot = intersect(c("S100A8", "S100A9", "IFITM2", "FCGR3B"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Basophil / Mast
FeaturePlot(object = expr, features.plot = intersect(c("MS4A2", "CPA3", "TPSAB1"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Eosinophil
FeaturePlot(object = expr, features.plot = intersect(c("SIGLEC8"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Megakaryocyte
FeaturePlot(object = expr, features.plot = intersect(c("NRGN", "PPBP", "PF4", "OST4"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Macrophage
FeaturePlot(object = expr, features.plot = intersect(c("MARCO", "MSR1", "MRC1"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# DC
FeaturePlot(object = expr, features.plot = intersect(c("HLA-DRB1", "CLEC9A", "LAMP3", "CD1C", "PLD4"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Monocyte
FeaturePlot(object = expr, features.plot = intersect(c("CD14", "S100A8", "FCGR3A"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

FeaturePlot(object = expr, features.plot = intersect(c("CD3D", "CD3E", "TOP2A", "MKI67"), rownames(expr@scale.data)), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

x <- "COL1A1"
do_queryGene(x)
FeaturePlot(object = expr, features.plot = x, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.hover = T)
VlnPlot(object = expr, features.plot = c(x), size.title.use = 14, point.size.use = 0.1)

xxx <- FindMarkers(object = expr, ident.1 = 0, ident.2 = NULL, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx)

yyy <- FindMarkers(object = expr, ident.1 = 2, ident.2 = NULL, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(yyy)

zzz <- FindMarkers(object = expr, ident.1 = 4, ident.2 = NULL, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(zzz)

kkk <- FindMarkers(object = expr, ident.1 = 12, ident.2 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(kkk)
kkk <- subset(kkk, power>=0.4 & avg_logFC>=log(2))
cat(rownames(kkk), sep = "\n")

kkk <- FindMarkers(object = expr, ident.1 = 10, ident.2 = 7, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
kkk <- subset(kkk, power>=0.4 & avg_logFC>=log(1.5))
cat(rownames(kkk), sep = "\n")

sss <- FindMarkers(object = expr, ident.1 = 3, ident.2 = 5, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss)

xxx <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, left1 = -4, right1 = -1, bottom1 = -18, top1 = -16.5)
head(xxx)

yyy <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = -11, right1 = -5, bottom1 = 9.5, top1 = 12.5)
head(yyy)

tmp <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_unfiltered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
dim(tmp)
tmp <- tmp[, colnames(expr_data)]
xxx <- tmp[grep("^MUC", rownames(tmp)), ]
