#!/usr/bin/env r
suppressMessages({
  library(docopt)
})


# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [--root ROOT]

-h --help           show help message
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(gridExtra)
  library(cicero)
  # library(networkD3)
  library(plyr)
  
  library(miniUI)
  library(shiny)
})

root <- opt$root
stage <- "11-14w"

setwd(paste(root, "data/11-14w/15_liver/results", sep = "/"))

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

proj <- loadArchRProject("ArchR/ArchR_output/", showLogo = F)

# remove cells in S5_plate3
# proj <- proj[!proj$plate %in% "PD10_HCA_15_12-14w_S5_Plate3"]


# de novo clustering
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
# print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))

proj <- addClusters(proj, reducedDims = "peakLSI", name = "clusters", dimsToUse = 1:30, force = T, k.param = 15, prune.SNN = 1/15, resolution = 0.3)

print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))

proj$tuned_ident <- mapvalues(proj$clusters, from = paste0("C", 1:11),
                              to = c("Erythrocyte", "Erythrocyte", "Erythrocyte", "Macro", "Hepatocyte", 
                                     "Endo", "Fibro", "T", "B", "Epi", 
                                     "Macro"))


print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_ident", plotAs = "points", size = 1.5))

# manual curation of clusters, two version, (1) fixed coordination; (2) shiny interactive session ----
# interactive_choose_cell(proj, c("EPCAM", "MARCO", "CD68", "CYP3A7"), group_bys = c("tuned_ident", "plate"), log2Norm = T)

un_cell_names <- static_choose_cell(proj, xlimit = c(-9.57, -8.69), ylimit = c(1.20, 3.04), group_by = "tuned_ident", used_groups = "Epi")
un_cell_names <- static_choose_cell(proj, xlimit = c(-8.95, -8.56), ylimit = c(0.59, 1.03), group_by = "tuned_ident", used_groups = "Epi") %>% c(un_cell_names)
un_cell_names <- static_choose_cell(proj, xlimit = c(-9.02, -8.65), ylimit = c(1.18, 1.50), group_by = "tuned_ident", used_groups = "Epi") %>% c(un_cell_names)

un_cell_names <- static_choose_cell(proj, xlimit = c(-9.93, -8.50), ylimit = c(3.18, 4.26), group_by = "tuned_ident", used_groups = "Macro") %>% c(un_cell_names) # MPP like


tuned_idents <- getCellColData(proj, select = "tuned_ident")
tuned_idents[un_cell_names, ] <- "Unknown"
proj <- addCellColData(proj, data = tuned_idents[[1]], name = "tuned_ident", cells = rownames(tuned_idents), force = T)

print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_ident", plotAs = "points", size = 1.5))


# DE analysis of unknown clusters ----------------------------------------------
geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix") %>% 
  summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")

matDR <- getReducedDims(
  ArchRProj = proj, 
  reducedDims = "peakLSI", 
  dimsToUse = 1:30, 
  corCutOff = 0.75, 
  scaleDims = NULL
)

geneScore[['lsi']] <- Seurat::CreateDimReducObject(embeddings = matDR, key = 'LSI_', assay = 'ACTIVITY')
geneScore <- FindVariableFeatures(geneScore, nfeatures = 3000)

Idents(geneScore) <- geneScore[["tuned_ident"]]

geneScore <- RunTSNE(geneScore, reduction = "lsi", dims = seq_len(Embeddings(geneScore, "lsi") %>% ncol))
geneScore <- RunUMAP(geneScore, reduction = "lsi", dims = seq_len(Embeddings(geneScore, "lsi") %>% ncol), min.dist = 0.5)

markers <- FindAllMarkers(geneScore, assay = "ACTIVITY", test.use = "roc", features = VariableFeatures(geneScore), only.pos = T)
# markers <- FindAllMarkers(geneScore, assay = "ACTIVITY", test.use = "roc", only.pos = T)

topn_df <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)

pdf("markers/tuned_signature_genes.pdf")
tunedSignatureGenesPlot(proj, topn_df, geneScore)
dev.off()

markers <- markers %>% filter(myAUC >= 0.5)
markers$gene <- standardize_genename(markers$gene, getGeneAnnotation(proj)$genes$symbol)
write.table(markers, file = "markers/tuned_signature_genes.txt", quote = F, 
            sep = "\t", row.names = F)

pdf("markers/MPP.pdf")
plotMarker(proj, c("CD68"), do_plot = F)
dev.off()

pdf("markers/hepatocyte.pdf")
plotMarker(proj, c("CYP3A7"), do_plot = F)
dev.off()

pdf("markers/neutrophil.pdf")
plotMarker(proj, c("S100A8", "S100A9", "IFITM2", "FCGR3B"), do_plot = F)
dev.off()


# save results
proj$tuned_group <- ident2clgrp(proj$tuned_ident)

cellMeta <- getCellColData(proj)[, c("seqID", "tissue", "samplingPos", "plate", "individual", "PassQC", "predictedIdent", "Reads", "Aligned_ratio",
                                     "nFrags", "FRIP", "mito_ratio", "BlacklistRatio", "DoubletScore", "nPeak", "group", "tuned_group")]
colnames(cellMeta) <- c("seqID", "tissue", "samplingPos", "plate", "individual", "QC", "ident", "cleanReads", "mpRatio",
                        "nFrags", "FRIP", "MitoRatio", "BlacklistRatio", "DoubletScore", "nPeak", "group", "tuned_group")

cellMeta$stage <- stage
cellMeta$species <- "human"
cellMeta$QC <- as.logical(cellMeta$QC)

# add umap/tsne embedding
peakTSNE <- getEmbedding(proj, embedding = "peakTSNE")
peakUMAP <- getEmbedding(proj, embedding = "peakUMAP")

colnames(peakTSNE) <- c("tSNE_1", "tSNE_2")
colnames(peakUMAP) <- c("UMAP_1", "UMAP_2")

cellMeta <- cbind(cellMeta, peakTSNE)
cellMeta <- cbind(cellMeta, peakUMAP)

write.table(cellMeta, file = "tuned_filtered_cellMeta_internal.txt", sep = "\t", quote = F, col.names = NA)

rownames(cellMeta) <- gsub("^.*#", "", rownames(cellMeta))

if (stage == "19-22w") {
  rownames(cellMeta) <- gsub("^(.*?_)", "\\1A_", rownames(cellMeta)) 
} else{
  rownames(cellMeta) <- gsub("^(.*?_)", "\\1B_", rownames(cellMeta))
}

write.table(cellMeta, file = "tuned_filtered_cellMeta.txt", sep = "\t", quote = F, col.names = NA)

saveArchRProject(proj, load = F)
