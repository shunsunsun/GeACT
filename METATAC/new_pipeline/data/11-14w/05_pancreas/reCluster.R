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
})

root <- opt$root
stage <- "11-14w"

setwd(paste(root, "data/11-14w/05_pancreas/results", sep = "/"))

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

proj <- loadArchRProject("ArchR/ArchR_output/", showLogo = F)

# de novo clustering
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))

proj <- addClusters(proj, reducedDims = "peakLSI", name = "clusters", dimsToUse = 1:30, force = T, k.param = 10, resolution = 0.15)

print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))

proj$tuned_ident <- mapvalues(proj$clusters, from = paste0("C", 1:11),
                              to = c("Macro", "B", "T", "Erythrocyte", "Fibro", 
                                     "Fibro-like", "DPP6+PDX1-", "Glial", "Endo", "DPP6+PDX1+", 
                                     "Epi"))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_ident", plotAs = "points", size = 1.5))

# DE analysis of unknown clusters
geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix") %>% 
  summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")

matDR <- getReducedDims(
  ArchRProj = proj,
  reducedDims = "peakLSI",
  dimsToUse = 1:30,
  corCutOff = 0.75,
  scaleDims = NULL
)

geneScore[['lsi']] <- Seurat::CreateDimReducObject(embeddings=matDR, key='LSI_', assay='ACTIVITY')
geneScore <- FindVariableFeatures(geneScore, nfeatures = 3000)

Idents(geneScore) <- geneScore[["tuned_ident"]]


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

pdf("markers/alpha.pdf")
plotMarker(proj, genes = c("GCG", "LOXL4", "PLCE1", "IRX2", "GC", "KLHL41"), do_plot = F)
dev.off()

pdf("markers/beta.pdf")
plotMarker(proj, genes = c("INS", "IAPP", "MAFA", "NPTX2", "DLK1"), do_plot = F)
dev.off()

pdf("markers/others.pdf")
plotMarker(proj, genes = c("DPP6", "PDX1"), do_plot = F)
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
