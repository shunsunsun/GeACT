
TSNEPlot(object = expr, pt.size = 1, do.label = T, no.legend = T)

xxx <- read.table("03-expression/merged/filtering/filtering_cells.txt", header = T, sep = "\t", stringsAsFactors = F)
ggplot(subset(xxx, filter), aes(x = nUMI, y = nGene, color = X %in% "XCX_C-A6_01")) + geom_point(alpha = 0.6)

FeaturePlot(object = expr, features.plot = c("EPCAM", "VIM", "PTPRC", "HBG1", "PECAM1", "COL1A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F)
# smooth muscle
FeaturePlot(object = expr, features.plot = c("CNN1", "ACTA2", "TAGLN", "NOTCH3", "LGR5", "DES", "LGR6", "MYH11"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# fibroblast
FeaturePlot(object = expr, features.plot = c("COL1A1", "PDGFRA", "ELN", "ACTA2", "PLIN2", "APOE"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# other
FeaturePlot(object = expr, features.plot = c("CSPG4", "TRPC6", "PDGFRB"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

## Podocyte is EPCAM negative https://www.karger.com/Article/FullText/111039
FeaturePlot(object = expr, features.plot = c("EPCAM"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
##

## Science (by Young)
# pRCC (papillary renal cell carcinoma)
#FeaturePlot(object = expr, features.plot = c("MET"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# ccRCC (clear cell renal cell carcinoma)
#FeaturePlot(object = expr, features.plot = c("CA9", "NDUFA4L2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# cap mesenchyme (CM)
FeaturePlot(object = expr, features.plot = c("SIX2", "CITED1", "PAX2", "SIX1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Collecting_duct (NA, NA, B, A)
FeaturePlot(object = expr, features.plot = c("ATP6V0D2", "CLCNKB", "SLC26A4", "SLC4A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Distal
FeaturePlot(object = expr, features.plot = c("AVPR2", "SLC8A1", "KCNJ1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Distal_collecting
FeaturePlot(object = expr, features.plot = c("CLDN8"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Glom_vascular
FeaturePlot(object = expr, features.plot = c("CLDN5", "SEMA3G", "AQP1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Glomerulus (Podocyte)
FeaturePlot(object = expr, features.plot = c("PTPRO", "PODXL", "WT1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Henle_ascending
FeaturePlot(object = expr, features.plot = c("CLDN16"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Henle_descending
FeaturePlot(object = expr, features.plot = c("SLC12A1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Myofibroblasts
FeaturePlot(object = expr, features.plot = c("PDGFRB", "ACTA2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Pelvic
FeaturePlot(object = expr, features.plot = c("KRT23", "SAA2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Pelvic_ureter (cluster 6)
FeaturePlot(object = expr, features.plot = c("TP63", "KRT5", "UPK1B", "UPK1A", "DHRS2", "S100P"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Principal_cells
FeaturePlot(object = expr, features.plot = c("AQP2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Proximal
FeaturePlot(object = expr, features.plot = c("SLC13A3", "SLC34A1", "SLC17A3", "SLC22A8", "SLC7A13", "SLC16A9", "SLC22A7"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# UretericBud
FeaturePlot(object = expr, features.plot = c("HNF1B", "RET", "GATA3", "ELF3", "POU3F3", "TFCP2L1", "CDH16"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# Vascular
FeaturePlot(object = expr, features.plot = c("PLVAP", "SLC14A1", "VCAM1", "KDR", "PTPRB", "PECAM1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

## Cell report
# CM
FeaturePlot(object = expr, features.plot = c("SIX2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# RI
FeaturePlot(object = expr, features.plot = c("SFRP1", "MEIS1", "PDGFRA"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# collecting duct
FeaturePlot(object = expr, features.plot = c("AQP2", "AQP3"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# PT
FeaturePlot(object = expr, features.plot = c("LRP2", "CUBN"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# podocyte
FeaturePlot(object = expr, features.plot = c("NPHS1", "NPHS2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# loop of Henle
FeaturePlot(object = expr, features.plot = c("UMOD", "POU3F3"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# distal convoluted tubule (DT)
FeaturePlot(object = expr, features.plot = c("CLCNKB"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# extraglomerular mesangium (EM)
FeaturePlot(object = expr, features.plot = c("PDGFRB", "REN"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
#  intraglomerular mesangium (MG)
FeaturePlot(object = expr, features.plot = c("PDGFRB", "LAMA4"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

## PNAS
# PC
FeaturePlot(object = expr, features.plot = c("AQP2", "AQP3", "SCNN1B", "SCNN1G", "KCNJ1", "AVPR2"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)
# IC
FeaturePlot(object = expr, features.plot = c("CA2", "SLC4A1", "SLC26A4", "RHCG", "ATP6V1B1"), cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = T)

x <- "S100B"
do_queryGene(x)
FeaturePlot(object = expr, features.plot = x, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2, no.legend = F, do.hover = T)
VlnPlot(object = expr, features.plot = c(x), size.title.use = 14, point.size.use = 0.1)

xxx1 <- FindMarkers(object = expr, ident.1 = 0, ident.2 = c(17,21), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx1)

xxx2 <- FindMarkers(object = expr, ident.1 = 17, ident.2 = c(0,21), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(xxx2)

yyy <- FindMarkers(object = expr, ident.1 = c(6,18,13,25), ident.2 = c(20,22,8,10), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(yyy)

zzz <- FindMarkers(object = expr, ident.1 = 0, ident.2 = c(16, 19), test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(zzz)

kkk <- FindMarkers(object = expr, ident.1 = 1, ident.2 = 5, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(kkk)

sss <- FindMarkers(object = expr, ident.1 = 5, ident.2 = 1, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss)

sss1 <- FindMarkers(object = expr, ident.1 = 11, ident.2 = 2, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss1)

sss2 <- FindMarkers(object = expr, ident.1 = 2, ident.2 = 11, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
head(sss2)

xxx <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = -14.2, right1 = -10, top1 = 16, bottom1 = 11, 
                       left2 = -10, right2= 5, top2 = -44, bottom2 = -58
                       )
head(xxx)

zzz <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = 18, right1 = 23, top1 = -7, bottom1 = -11)
head(zzz)

kkk <- do_regionMarker(expr = expr, cellMetaData = cellMetaData, 
                       left1 = -5, right1 = 5, bottom1 = -36, top1 = -33, 
                       left2 = -5, right2 = 5, bottom2 = -33, top2 = -20)
head(kkk)

tmp <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_unfiltered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
dim(tmp)
tmp <- tmp[, colnames(expr_data)]
xxx <- tmp[grep("^MUC", rownames(tmp)), ]
