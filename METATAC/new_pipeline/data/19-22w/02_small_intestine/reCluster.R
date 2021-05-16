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
setwd(paste(root, "data/19-22w/02_small_intestine/results", sep = "/"))

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

proj <- loadArchRProject("ArchR/ArchR_output/", showLogo = F)

# proj <- addUMAP(
#   ArchRProj = proj,
#   reducedDims = "peakLSI",
#   name = "peakUMAP",
#   nNeighbors = 30,
#   minDist = 0.6,
#   metric = "cosine",
#   force = T
# )



print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
# print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))

proj <- addClusters(proj, reducedDims = "peakLSI", name = "clusters", dimsToUse = 1:30, force = T, k.param = 15, prune.SNN = 1/15, resolution = 0.3)

print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakTSNE", colorBy = "cellColData", name = "clusters", plotAs = "points", size = 1.5))

proj$tuned_label <- mapvalues(proj$clusters, from = paste0("C", 1:12),
                              to = c("Macrophage", "T", "B", "Fibroblast", "Fibroblast", 
                                     "Smooth muscle", "Fibroblast", "Fibroblast", "Fibroblast", 
                                     "Epithelial", "Glial", "Endothelial"))


print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_label", plotAs = "points", size = 1.5))

# manual curation of clusters, two version, (1) fixed coordination; (2) shiny interactive session ----
interactive_choose_cell(proj, "EPCAM", log2Norm = T, group_bys = "tuned_label")

un_cell_names  <- static_choose_cell(proj, xlimit = c(-0.21, 3.93), ylimit = c(0.92, 5.07), group_by = "tuned_label", used_groups = "Epithelial")

tuned_labels <- getCellColData(proj, select = "tuned_label")
tuned_labels[un_cell_names, ] <- "Unknown"

proj <- addCellColData(proj, data = tuned_labels[[1]], name = "tuned_label", cells = rownames(tuned_labels), force = T)

print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_label", plotAs = "points", size = 1.5))


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

Idents(geneScore) <- geneScore[["tuned_label"]]

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


# helper functions, not run ----------------------------------------------------
# confusion matrix to exclude batch clusters
cM <- confusionMatrix(paste0(proj$clusters), paste0(proj$plate))
cM
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p


# SNN graph on tSNE to explain some weird clustering results
matDR <- getReducedDims(
  ArchRProj = proj, 
  reducedDims = "peakLSI", 
  dimsToUse = 1:30, 
  corCutOff = 0.75, 
  scaleDims = NULL
)

tmp <- matrix(rnorm(nrow(matDR) * 3, 10), ncol = nrow(matDR), nrow = 3)
colnames(tmp) <- rownames(matDR)
rownames(tmp) <- paste0("t", seq_len(nrow(tmp)))

obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
obj[['lsi']] <- Seurat::CreateDimReducObject(embeddings=matDR, key='LSI_', assay='RNA')

obj <- RunTSNE(obj, reduction = "lsi", dims = seq_len(ncol(matDR)))
obj <- RunUMAP(obj, reduction = "lsi", dims = seq_len(ncol(matDR)), min.dist = 0.5)

obj <- FindNeighbors(obj, dims = seq_len(ncol(matDR)), reduction = "lsi", do.plot = T, k.param = 15, prune.SNN = 1/5)

obj <- FindClusters(obj, resolution = 0.1)

DimPlot(obj, label = T, reduction = "tsne")
DimPlot(obj, label = T)

proj$clusters <- Idents(obj)[getCellNames(proj)]




interactive_choose_df(data_df, feature = "HBB", point_size = 1.5, group_bys = c("plate", "tuned_label"))
  
library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

