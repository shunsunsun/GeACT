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

groups_to_keep <- c("Epithelial", "Endothelial", "Fibroblast", "Glial", "T")
tissue_order <- c("esophagus", "stomach", "small intestine", "large intestine", "pancreas", "liver", "kidney", "bladder", "lung", "bronchus", "diaphragm", 
                  "spleen", "bone marrow", "thymus", "heart", "ovary", "testis")

dir.create("fig5/supp", showWarnings = F, recursive = T)

for (group in groups_to_keep) {
  if( dir.exists(paste0("fig5/supp/ArchR_", group)) ){
    proj_group <- loadArchRProject(paste0("fig5/supp/ArchR_", group), showLogo = F)
  } else {
    cat("Calculating marker motifs for ", group, "\n")
    
    proj_group <- proj_19_22w[proj_19_22w$tuned_group == group, ]
    proj_group <- saveArchRProject(proj_group, outputDirectory = paste0("fig5/supp/ArchR_", group), load = T, dropCells = T)
    
    proj_group <- addMotifAnnotations(ArchRProj = proj_group, motifSet = use_motif_set, name = paste0(use_motif_set, "Motif"), force = T)
    proj_group <- addBgdPeaks(proj_group, method = "ArchR", force = T)
    proj_group <- addDeviationsMatrix(
      ArchRProj = proj_group, 
      peakAnnotation = paste0(use_motif_set, "Motif"),
      force = TRUE
      )
    
    proj_group <- saveArchRProject(proj_group, load = T)
    
    motif_matrix <- getMatrixFromProject(proj_group, useMatrix = paste0(use_motif_set, "MotifMatrix"))
    rownames(motif_matrix) <- gsub("_.*", "", rownames(motif_matrix))
    
    # # standardize JASPAR motif names
    motif_matrix_zscore <- assays(motif_matrix)$z
    if (use_motif_set == "JASPAR2020") {
      rownames(motif_matrix_zscore) <- rownames(motif_matrix_zscore) %>% gsub("\\.\\.", "::", .) %>% 
        gsub("\\.(var\\.[0-9]*)$", "(\\1)", .) %>%
        gsub("(?<!var)\\.", "-", ., perl = T)
    }
    
    motif_matrix.seurat <- CreateSeuratObject(counts = motif_matrix_zscore, project = "chromVAR", assay = "TF", names.delim = "#", meta.data = as.data.frame(colData(motif_matrix)))
    motif_matrix.seurat <- ScaleData(motif_matrix.seurat)
    Idents(motif_matrix.seurat) <- "tissue"
    
    Idents(motif_matrix.seurat) <- first_case_up(Idents(motif_matrix.seurat) %>% as.character)
    levels(motif_matrix.seurat) <- levels(motif_matrix.seurat)[match(levels(motif_matrix.seurat), first_case_up(tissue_order) ) %>% order]
    
    motif_matrix.seurat <- RunPCA(motif_matrix.seurat, features = rownames(motif_matrix.seurat))
    motif_matrix.seurat <- RunUMAP(motif_matrix.seurat, dims = 1:30)
    DimPlot(motif_matrix.seurat, label = T)
    
    marker_motifs <- FindAllMarkers(motif_matrix.seurat, test.use = "roc", only.pos = TRUE, min.pct = 0.1)
    marker_motifs$motif <- gsub("-.*", "", marker_motifs$gene)
    
    # change this remove process to group level
    marker_motifs <- marker_motifs %>% group_by(cluster)
    
    write.table(marker_motifs, file = paste0("fig5/supp/", group, "_marker_motifs.txt"), sep = "\t", quote = F, row.names = F)
    
    marker_motifs <- group_map(marker_motifs, ~.x[!duplicated(.x$motif), ], .keep = T) %>% bind_rows()
    
    write(paste0(group, ": ", length(unique(marker_motifs$motif))), "fig5/supp/summary.txt", append = T)
  }
}
