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
stage <- "19-22w"

setwd(paste(root, "data/19-22w/01_stomach/results", sep = "/"))

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

proj <- loadArchRProject("ArchR/ArchR_output/", showLogo = F)


# de novo clustering
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
# print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))

proj <- addClusters(proj, reducedDims = "peakLSI", name = "clusters", dimsToUse = 1:30, force = T, k.param = 10, prune.SNN = 1/15, resolution = 0.15)

print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))

proj$tuned_ident <- mapvalues(proj$clusters, from = paste0("C", 1:11),
                              to = c("Fibro", "Fibro", "Fibro", "Fibro", "Fibro", 
                                     "Glial", "Unknown", "Epi", "Endo", "T", 
                                     "Macro"))


print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_ident", plotAs = "points", size = 1.5))

# manual curation of clusters, two version, (1) fixed coordination; (2) shiny interactive session ----
# interactive_choose_cell(proj, c("PECAM1", "MS4A1", "MARCO", "HBG2", "CACNA1A", "CD3D"), log2Norm = T, group_bys = c("tuned_ident", "plate"))


un_cell_names  <- c(static_choose_cell(proj, xlimit = c(-9.06, -8.47), ylimit = c(10.67, 11.09)), 
                    static_choose_cell(proj, xlimit = c(2.16, 2.42), ylimit = c(11.05, 11.43)))
B_cell_names   <- static_choose_cell(proj, xlimit = c(2.37, 2.97), ylimit = c(10.97, 11.76))
ery_cell_names <- static_choose_cell(proj, xlimit = c(1.73, 3.27), ylimit = c(11.80, 13.52))


tuned_idents <- getCellColData(proj, select = "tuned_ident")
tuned_idents[un_cell_names, ] <- "Unknown"
tuned_idents[B_cell_names, ] <- "B"
tuned_idents[ery_cell_names, ] <- "Erythrocyte"

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

topn_df$gene <- standardize_genename(topn_df$gene, getGeneAnnotation(proj)$genes$symbol)
write.table(topn_df, file = "markers/tuned_signature_genes.txt", quote = F, 
            sep = "\t", row.names = F)

markers <- markers %>% filter(myAUC >= 0.5)
markers$gene <- standardize_genename(markers$gene, getGeneAnnotation(proj)$genes$symbol)
write.table(markers, file = "markers/tuned_signature_genes.txt", quote = F, 
            sep = "\t", row.names = F)

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
