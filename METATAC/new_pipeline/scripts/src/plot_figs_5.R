#!/usr/bin/env r
suppressPackageStartupMessages({
  library(docopt)
})

# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [--root ROOT] [-m MOTIFSET]

-h --help           show help message
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
-m --motif MOTIFSET motif set to use [default: JASPAR2020]
"

opt <- docopt(doc)

root <- opt$root
use_motif_set <- opt$motif

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
  library(circlize)
  library(GenomicRanges)
  library(RColorBrewer)
  library(ComplexHeatmap)
})

setwd("/data/Lab/otherwork/GeACT/ATAC/data/all")
dir.create("fig5", showWarnings = F)
set.seed(seed = 0)

proj <- loadArchRProject(path = "ArchR_output/", showLogo = F)

if(dir.exists("ArchR_19_22w")){
  proj_19_22w <- loadArchRProject("ArchR_19_22w", showLogo = F)
} else{
  proj_19_22w <- proj[proj$stage == "19-22w", ]
  
  proj_19_22w <- saveArchRProject(proj_19_22w, outputDirectory = "ArchR_19_22w", load = T, dropCells = T)
  
  proj_19_22w <- addIterativeLSI(proj_19_22w, useMatrix = "PeakMatrix", 
                              iterations = 1, name = "peakLSI", varFeatures = 50000, force = T)
  
  proj_19_22w <- addTSNE(
    ArchRProj = proj_19_22w, 
    reducedDims = "peakLSI", 
    name = "peakTSNE", 
    perplexity = 30, 
    force = T
  )
  
  proj_19_22w <- addUMAP(
    ArchRProj = proj_19_22w, 
    reducedDims = "peakLSI", 
    name = "peakUMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", 
    force = T
  )
}

# fig 5 ------------------------------------------------------------------------

# 1 Compare same group of cells across different tissues -----------------------
if(dir.exists("ArchR_epi")){
  proj_epi <- loadArchRProject("ArchR_epi", showLogo = F)
} else {
  proj_epi <- proj_19_22w[proj_19_22w$tuned_group == "Epithelial", ]
  
  proj_epi <- saveArchRProject(proj_epi, outputDirectory = "ArchR_epi", load = T, dropCells = T)
  
  # 1.1 visualize clusters -----------------------------------------------------
  proj_epi <- addIterativeLSI(proj_epi, useMatrix = "PeakMatrix", 
                              iterations = 1, name = "peakLSI", varFeatures = 50000, force = T)
  
  proj_epi <- addTSNE(
    ArchRProj = proj_epi, 
    reducedDims = "peakLSI", 
    name = "peakTSNE", 
    perplexity = 30, 
    force = T
  )
  
  proj_epi <- addUMAP(
    ArchRProj = proj_epi, 
    reducedDims = "peakLSI", 
    name = "peakUMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", 
    force = T
  )
  
  p1 <- plotEmbedding(proj_epi, embedding = "peakUMAP", colorBy = "cellColData", name = "tissue", size = 1)
  pdf("fig5/dimred.pdf", width = 10, height = 10)
  p1
  dev.off()
  
  cell_meta <- getCellColData(proj_epi, select = "tissue")
  dimred <- getEmbedding(proj_epi, embedding = "peakUMAP")
  
  df_umap <- merge(cell_meta, dimred, by = 0) %>% as.data.frame %>% column_to_rownames("Row.names")
  names(df_umap) <- c("tissue", "UMAP_1", "UMAP_2")
  
  write.table(df_umap, file = "epi_dimred.txt", sep = "\t", quote = F, col.names = NA)
}



# 1.2 visualize peaks and peak-gene linkage aside marker genes -----------------

# # coaccessibility
# proj_epi <- addCoAccessibility(proj_epi, reducedDims = "peakLSI", dimsToUse = 1:30)
# 
# # peak to gene linkage
# rna_files <- dir(paste(root, "data/19-22w", sep = "/"), pattern = "expr_RNA.rds", recursive = T, full.names = T)
# rna_files <- grep("RNA/19-22w/expr_RNA.rds", rna_files, value = T)
# # rna_files <- grep("06_spleen|07_testis", rna_files, value = T, invert = T)
# tissues <- gsub(paste(root, "data/19-22w/", sep = "/"), "", rna_files) %>% gsub("/.*", "", .) %>% gsub("^.*?_", "", .) %>% gsub("_", " ", .)
# names(rna_files) <- tissues
# 
# my_merge <- function(s1, s2){
#   if (is.null(s1))
#     return(s2)
#   if (is.null(s2))
#     return(s1)
#   return(merge(s1, s2))
# }
# 
# epi_seurat_obj <- lapply(seq_along(tissues), function(i) {
#   sobj <- readRDS(rna_files[i])
#   sobj$group <- ident2clgrp(sobj$ident)
#   sobj$tissue <- tissues[i]
#   cat(tissues[i], "\n")
#   print(unique(sobj$group))
#   if ("Epithelial" %in% unique(sobj$group)) {
#     return(sobj[, sobj$group == "Epithelial"])
#   } else {
#     return(NULL)
#   }
# }) %>% reduce(my_merge)
# 
# epi_seurat_obj <- epi_seurat_obj[, epi_seurat_obj$tissue %in% unique(proj_epi$tissue)]
# 
# proj_epi <- addGeneScoreMatrix(proj_epi, force = T)
# 
# proj_epi <- addGeneIntegrationMatrix(proj_epi, useMatrix = "GeneScoreMatrix",
#                                      reducedDims = "peakLSI", seRNA = epi_seurat_obj,
#                                      addToArrow = TRUE,
#                                      force= TRUE,
#                                      groupRNA = "tissue",
#                                      nameCell = "predictedCell",
#                                      nameGroup = "predictedTissue",
#                                      nameScore = "predictedScore", 
#                                      plotUMAP = T, 
#                                      useImputation = F)
# 
# proj_epi <- addPeak2GeneLinks(proj_epi,
#                               reducedDims = "peakLSI")


# genome track for marker genes in epi of different tissues 
markerGenes  <- c(
  "KRT15", # Esophagus
  "CLDN18", # Stomach
  "MUC13", # Small intestine
  "LGALS4", # Large intestine
  "CLPS", # Pancreas
  "MYB", # Liver
  "PODXL", # Kidney
  "CPM" # Lung
  # "ANXA1", # Bronchus, Esophagus
  # "CCDC170" # Ovary
)

p <- plotBrowserTrack(
  ArchRProj = proj_epi, 
  groupBy = "tissue", 
  geneSymbol = markerGenes, 
  upstream = 100000,
  downstream = 100000,
  baseSize = 12,
  facetbaseSize = 12,
  loops = NULL # getCoAccessibility(proj_epi) # getPeak2GeneLinks(proj_epi)
)

pdf("fig5/epi_tissue_gene_markers_genome_track.pdf", width = 10, height = 10)
for (fig in p){
  grid::grid.newpage()
  grid::grid.draw(fig)
}
dev.off()





# 1.3 chromVAR TF motif analysis -----------------------------------------------
proj_epi <- addMotifAnnotations(ArchRProj = proj_epi, motifSet = use_motif_set, name = paste0(use_motif_set, "Motif"), force = T)
proj_epi <- addBgdPeaks(proj_epi, method = "ArchR")
proj_epi <- addDeviationsMatrix(
  ArchRProj = proj_epi, 
  peakAnnotation = paste0(use_motif_set, "Motif"),
  force = TRUE
)

proj_epi <- saveArchRProject(proj_epi, load = T)
# proj_epi <- loadArchRProject("ArchR_epi")


plotVarDev <- getVarDeviations(proj_epi, name = paste0(use_motif_set, "MotifMatrix"), plot = TRUE)

pdf(paste0("fig5/epi_tissue_", use_motif_set, "_motif_markers.pdf"), width = 10, height = 10)
plotVarDev

# 1.3.1 Seurat analysis
motif_matrix <- getMatrixFromProject(proj_epi, useMatrix = paste0(use_motif_set, "MotifMatrix"))
rownames(motif_matrix) <- gsub("_.*", "", rownames(motif_matrix))

# motif_matrix.seurat <- as.Seurat(motif_matrix, counts = "deviations", data = "z", project = "chromVAR")
# colnames(motif_matrix)
motif_matrix.seurat <- CreateSeuratObject(counts = assays(motif_matrix)$z, project = "chromVAR", assay = "TF", names.delim = "#", meta.data = as.data.frame(colData(motif_matrix)))
motif_matrix.seurat <- ScaleData(motif_matrix.seurat)
Idents(motif_matrix.seurat) <- "tissue"
levels(motif_matrix.seurat) <- c("esophagus", "stomach", "small intestine", "large intestine", "pancreas", "liver", "kidney", "lung")

motif_matrix.seurat <- RunPCA(motif_matrix.seurat, features = rownames(motif_matrix.seurat))
motif_matrix.seurat <- RunUMAP(motif_matrix.seurat, dims = 1:30)
DimPlot(motif_matrix.seurat, label = T)

marker_motifs <- FindAllMarkers(motif_matrix.seurat, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
marker_motifs$motif <- gsub("-.*", "", marker_motifs$gene)
marker_motifs <- marker_motifs[!duplicated(marker_motifs$motif), ]

top8 <- marker_motifs %>% group_by(cluster) %>% top_n(n = 8, wt = myAUC) # p_val_adj)

DoHeatmap(motif_matrix.seurat, features = top8$gene, cells = WhichCells(motif_matrix.seurat, downsample = 200), angle = 60) + NoLegend() + theme(plot.margin = margin(t = 30, r = 40))
DoHeatmap(motif_matrix.seurat, features = top8$gene, angle = 60) + NoLegend() + theme(plot.margin = margin(t = 30, r = 40))
umap_epi <- getEmbedding(proj_epi, embedding = "peakUMAP")
if (use_motif_set == "cisbp"){
  # lung
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "NKX24", data = umap_epi, size = 3))
  # esophagus
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "GSC", data = umap_epi, size = 3))
  # stomach
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "FOXC2", data = umap_epi, size = 3))
  # small intestine
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "CDX1", data = umap_epi, size = 3))
  # large intestine
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "HOXD13", data = umap_epi, size = 3))
  # liver
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "GATA4", data = umap_epi, size = 3))
  # pancreas
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "HOXB7", data = umap_epi, size = 3))
  # kidney
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "PAX1", data = umap_epi, size = 3))
} else if (use_motif_set == "JASPAR2020"){
  # lung
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "NKX2.8", data = umap_epi, size = 3))
  # esophagus
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "SOX4", data = umap_epi, size = 3))
  # stomach
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "FOXC1", data = umap_epi, size = 3))
  # small intestine
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "CDX1", data = umap_epi, size = 3))
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "CDX2", data = umap_epi, size = 3))
  # large intestine
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "HOXC13", data = umap_epi, size = 3))
  # liver
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "GATA1..TAL1", data = umap_epi, size = 3))
  # pancreas
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "ONECUT2", data = umap_epi, size = 3))
  # kidney
  print(chromVARFeaturePlot(object = motif_matrix.seurat, feature = "PAX1", data = umap_epi, size = 3))
}




# 1.3.2 ArchR analysis (!with bug)
# 2 questions: 
# why producing many NAs (logFC): because tf motif score contain minus values, which cause problems when calculating logs
# how seurat tackle with this problem?
# why TF order is changed during plotting: due to row cluster
marker_motifs <- getMarkerFeatures(
  ArchRProj = proj_epi, 
  useMatrix = paste0(use_motif_set, "MotifMatrix"), 
  groupBy = "tissue",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z"
)

rownames(marker_motifs) <- gsub("_.*", "", rowData(marker_motifs)$name)

marker_list <- getMarkers(marker_motifs, cutOff = "FDR <= 0.01 & Log2FC >= 1")
marker_list

marker_list2 <- getMarkers(marker_motifs, cutOff = "FDR <= 0.01")
marker_list2

idx <- lapply(marker_list2, function(df){
  return(head(df$idx, 5))
}) %>% unlist
names(idx) <- NULL
idx <- unique(idx)

marker_motifs_toplot <- marker_motifs[idx, ]

heatmap_motifs <- plotMarkerHeatmap(
  seMarker = marker_motifs_toplot, 
  cutOff = "FDR <= 0.01",
  log2Norm = F,
  transpose = F,
  binaryClusterRows = F,
  clusterCols = T
)

draw(heatmap_motifs, heatmap_legend_side = "bot", annotation_legend_side = "bot")

dev.off()



# 2 compare types of cell from all tissues -------------------------------------

# 2.1 heatmap ------------------------------------------------------------------

# load color palette
load(paste(root, "meta/cg_color.RData", sep = "/"))

proj_19_22w <- addMotifAnnotations(proj_19_22w, motifSet = use_motif_set, name = paste0(use_motif_set, "Motif"))
proj_19_22w <- addBgdPeaks(proj_19_22w, method = "ArchR")
proj_19_22w <- addDeviationsMatrix(
  ArchRProj = proj_19_22w, 
  peakAnnotation = paste0(use_motif_set, "Motif"),
  force = TRUE
)


proj_19_22w <- saveArchRProject(proj_19_22w, outputDirectory = "ArchR_19_22w", load = T)
# proj_19_22w <- loadArchRProject(path = "ArchR_19_22w", showLogo = F)

plotVarDev <- getVarDeviations(proj_19_22w, name = paste0(use_motif_set, "MotifMatrix"), plot = TRUE)

pdf(sprintf("fig5/group_%s_motif_markers.pdf", use_motif_set), width = 10, height = 10)
plotVarDev

# Seurat analysis
motif_matrix <- getMatrixFromProject(proj_19_22w, useMatrix = paste0(use_motif_set, "MotifMatrix"))
rownames(motif_matrix) <- gsub("_.*", "", rownames(motif_matrix))

motif_matrix.seurat <- CreateSeuratObject(counts = assays(motif_matrix)$z, project = "chromVAR", assay = "TF", names.delim = "#", meta.data = as.data.frame(colData(motif_matrix)))
motif_matrix.seurat <- ScaleData(motif_matrix.seurat)
Idents(motif_matrix.seurat) <- "group"

marker_motifs <- FindAllMarkers(motif_matrix.seurat, only.pos = TRUE, test.use = "roc", min.pct = 0.25)

groups_to_keep <- c("Epithelial", "Endothelial", "Fibroblast", "Glial", "T")   # "Smooth muscle" "Erythrocyte"
motif_matrix.seurat <- motif_matrix.seurat[, Idents(motif_matrix.seurat) %in% groups_to_keep]
levels(motif_matrix.seurat) <- groups_to_keep

marker_motifs <- marker_motifs[marker_motifs$cluster %in% groups_to_keep, ]
marker_motifs$cluster <- factor(marker_motifs$cluster, levels = groups_to_keep)
marker_motifs <- marker_motifs[order(marker_motifs$cluster), ]

top10 <- marker_motifs %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC) # p_val_adj)
DoHeatmap(motif_matrix.seurat, features = top10$gene, cells = WhichCells(motif_matrix.seurat, downsample = 200), group.colors = cg_color, angle = 60, raster = T) + NoLegend() +
  theme(plot.margin = margin(t = 20, r = 40))
DoHeatmap(motif_matrix.seurat, features = top10$gene, group.colors = cg_color, angle = 60, raster = T) + NoLegend() +
  theme(plot.margin = margin(t = 20, r = 40))
dev.off()


# experimental
pdf("fig5/experiment.pdf", width = 10, height = 10)
umap_all <- getEmbedding(proj_19_22w[proj_19_22w$group %in% groups_to_keep, ], embedding = "peakUMAP")
plotEmbedding(proj_19_22w[proj_19_22w$group %in% groups_to_keep, ], embedding = "peakUMAP", colorBy = "cellColData", name = "group")
# Epithelial
chromVARFeaturePlot(object = motif_matrix.seurat, feature = "ZEB1", data = umap_all, size = 3)
# Endothelial
chromVARFeaturePlot(object = motif_matrix.seurat, feature = "SOX8", data = umap_all, size = 3)
# Fibroblast
chromVARFeaturePlot(object = motif_matrix.seurat, feature = "TAL1..TCF3", data = umap_all, size = 3)
# Glial
chromVARFeaturePlot(object = motif_matrix.seurat, feature = "SOX4", data = umap_all, size = 3)
# Erythrocyte
chromVARFeaturePlot(object = motif_matrix.seurat, feature = "SREBF1.var.2", data = umap_all, size = 3)



p <- plotBrowserTrack(
  ArchRProj = proj_19_22w[proj_19_22w$group %in% groups_to_keep], 
  groupBy = "group", 
  geneSymbol = "EPCAM", 
  upstream = 50000,
  downstream = 50000,
  baseSize = 12,
  facetbaseSize = 12,
  loops = NULL # getCoAccessibility(proj_epi) # getPeak2GeneLinks(proj_epi)
)

for (fig in p){
  grid::grid.newpage()
  grid::grid.draw(fig)
}

proj_19_22w <- addGeneScoreMatrix(proj_19_22w, force = T)
plotEmbedding(proj_19_22w[proj_19_22w$group %in% groups_to_keep, ], embedding = "peakUMAP", colorBy = "GeneScoreMatrix", name = "HBG1", plotAs = "points", size = 1, continuousSet = "whiteBlue", imputeWeights = NULL)
plotEmbedding(proj_19_22w[proj_19_22w$group %in% groups_to_keep, ], embedding = "peakUMAP", colorBy = "cellColData", name = "group", imputeWeights = NULL)

dev.off()

# 2.2 circos plot --------------------------------------------------------------
peak_matrix <- getMatrixFromProject(proj_19_22w, useMatrix = "PeakMatrix")

groups_to_cmp <- groups_to_keep

# first calculate marker peaks, later visualization is restricted to marker peaks
markers_peaks <- getMarkerFeatures(
  ArchRProj = proj_19_22w,
  useMatrix = "PeakMatrix",
  groupBy = "group",
  useGroups = groups_to_cmp,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
cut_off <- "FDR <= 0.01 & Log2FC >= 1"
marker_peak_list <- getMarkers(markers_peaks, cutOff = cut_off)

for(i in names(marker_peak_list)) {
  marker_peak_list[[i]]$cluster <- i
}

marker_peak_list_combined <- rbind(unlist(marker_peak_list))
marker_peak_list_combined$idx <- NULL
marker_peak_list_combined$end <- marker_peak_list_combined$end - 1 # bed transformation
write.table(marker_peak_list_combined, file = "fig5/group_marker_peaks.bed", col.names = T, quote = F, row.names = F, sep = "\t")

# determine if marker peaks are specific to one group
# for (i in seq_along(groups_to_cmp)) {
#   marker_peak_list[[i]]$group <- groups_to_cmp[i]
#   group_spec <- rep(TRUE, nrow(marker_peak_list[[i]]))
#   for (j in seq_along(groups_to_cmp)) {
#     if (i == j) next
#     cat(sprintf("Compare %s with %s\n", groups_to_cmp[i], groups_to_cmp[j]))
#     group_spec_markers <- getMarkerFeatures(
#       ArchRProj = proj_19_22w[proj_19_22w$group %in% groups_to_cmp[c(i, j)]],
#       useMatrix = "PeakMatrix",
#       groupBy = "group",
#       useGroups = groups_to_cmp[c(i, j)],
#       bias = c("TSSEnrichment", "log10(nFrags)"),
#       testMethod = "wilcoxon"
#     )
#     group_spec_markers_list <- getMarkers(group_spec_markers, cutOff = cut_off)
#     
#     group_spec <- group_spec & (rownames(marker_peak_list[[i]]) %in% 
#                                   rownames(group_spec_markers_list[[1]]))
#   }
#   marker_peak_list[[i]]$spec <- group_spec
# }

assayNames <- names(assays(markers_peaks))
for (an in assayNames) {
  eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(markers_peaks)[['", an, "']]")))
}
passMat <- eval(parse(text="FDR <= 1e-2 & Log2FC >= 1"))
for(an in assayNames){
  eval(parse(text=paste0("rm(",an,")")))
}
idx <- which(rowSums(passMat) > 0)

peaks <- getPeakSet(proj_19_22w)
peaks_sub <- peaks[idx, ]

# write marker peaks
# peaksDF <- data.frame(chr=seqnames(peaks_sub), start=start(peaks_sub) - 1, end=end(peaks_sub))
# peaksDF <- cbind(peaksDF, mcols(peaks_sub))
# write.table(peaksDF, file = "marker_peaks.bed", col.names = T, quote = F, row.names = F, sep = "\t")

peak_matrix_sub <- subsetByOverlaps(peak_matrix, peaks_sub)

use_colors <- cg_color[groups_to_cmp]

pdf(file = "fig5/circos.pdf", width = 10, height = 10)
plotCircosFromRangeSE(peak_matrix_sub, groups_to_cmp = groups_to_cmp, use_colors = use_colors)
dev.off()

pdf(file = "fig5/circos2.pdf", width = 10, height = 10)
plotCircosFromRangeSE(peak_matrix_sub, groups_to_cmp = groups_to_cmp, use_colors = use_colors, chrs = 2, window.size = 5e4, max.clip = 10)
dev.off()

pdf(file = "fig5/circos3.pdf", width = 10, height = 10)
epcam_df <- data.frame(chr = "chr2", start = 47345158, end = 47387601, name = "EPCAM")
cytoband_epcam <- list(df = data.frame(V1 = "chr2",
                                       V2 = 47000000,
                                       V3 = 48000000, 
                                       V4 = "p33", 
                                       V5 = "gpos50"), 
                       chromosome = "chr1", 
                       chr.len = setNames("chr1", 249250621))
plotCircosFromRangeSE(peak_matrix_sub, groups_to_cmp = groups_to_cmp, use_colors = use_colors, chrs = 2, window.size = 1e3, max.clip = 3,
                      cytoband = cytoband_epcam, genes_df = epcam_df)
dev.off()

pdf(file = "fig5/circos4.pdf", width = 10, height = 10)
plotCircosFromRangeSE(peak_matrix, groups_to_cmp = groups_to_cmp, use_colors = use_colors, chrs = 2, window.size = 1e3, max.clip = 3,
                      cytoband = cytoband_epcam, genes_df = epcam_df)
dev.off()



# tmp utils --------------------------------------------------------------------
# my_addPeak2GeneLinks <- edit(addPeak2GeneLinks)
# environment(my_addPeak2GeneLinks) <- asNamespace('ArchR')
