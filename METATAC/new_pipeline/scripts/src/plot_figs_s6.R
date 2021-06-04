#!/usr/bin/env r
suppressPackageStartupMessages({
  library(docopt)
})

# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [--root ROOT]

-h --help           show help message
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
"

opt <- docopt(doc)

root <- opt$root

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
set.seed(seed = 0)



if(dir.exists("ArchR_epi")){
  proj_epi <- loadArchRProject("ArchR_epi", showLogo = F)
} else {
  proj_19_22w <- loadArchRProject(path = "ArchR_19_22w/", showLogo = F)
  proj_epi <- proj_19_22w[proj_19_22w$tuned_group == "Epithelial", ]
  
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
  pdf("dimred.pdf", width = 10, height = 10)
  p1
  dev.off()
  
  cell_meta <- getCellColData(proj_epi, select = "tissue")
  dimred <- getEmbedding(proj_epi, embedding = "peakUMAP")
  
  df_umap <- merge(cell_meta, dimred, by = 0) %>% as.data.frame %>% column_to_rownames("Row.names")
  names(df_umap) <- c("tissue", "UMAP_1", "UMAP_2")
  
  write.table(df_umap, file = "epi_dimred.txt", sep = "\t", quote = F, col.names = NA)
}

# peak to gene linkage
rna_files <- dir(paste(root, "data/19-22w", sep = "/"), pattern = "expr_RNA.rds", recursive = T, full.names = T)
rna_files <- grep("RNA/19-22w/expr_RNA.rds", rna_files, value = T)
# rna_files <- grep("06_spleen|07_testis", rna_files, value = T, invert = T)
tissues <- gsub(paste(root, "data/19-22w/", sep = "/"), "", rna_files) %>% gsub("/.*", "", .) %>% gsub("^.*?_", "", .) %>% gsub("_", " ", .)
names(rna_files) <- tissues

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

# figure 6 ---------------------------------------------------------------------

# peak to gene links in epithelial/fibroblast
proj_epi



# tissue_order <- c("esophagus", "stomach", "small intestine", "large intestine", "pancreas", "liver", "kidney", "bladder", "lung", "bronchus", "diaphragm", 
#                   "spleen", "bone marrow", "thymus", "heart", "ovary", "testis")
# tissues <- proj_epi$tissue %>% unique
# tissues <- tissues[match(tissues, tissue_order) %>% order]
tissues <- c("esophagus", "stomach", "small intestine", "large intestine", "pancreas", "liver", "kidney", "lung")
names(tissues) <- tissues

# preprocess of rna seurat objects for integration
epi_seurat_objs <- lapply(tissues, function(tissue) {
  sobj <- readRDS(rna_files[tissue])
  sobj$group <- ident2clgrp(sobj$ident)
  sobj$tissue <- tissue
  cat(tissue, "\n")
  
  sobj[, sobj$group == "Epithelial"]
})

# ps <- plot_peak2gene_link(proj_epi, genes = markerGenes, outdir = "peak2gene_genescore", use_groups = tissues, groupby = "tissue", 
#                          use_matrix = "GeneScoreMatrix", save_p2g_mat = T, gene_size = 3)
ps <- plot_peak2gene_link(proj_epi, genes = markerGenes, outdir = "peak2gene_gs", use_groups = tissues, groupby = "tissue", 
                          use_matrix = "GeneScoreMatrix", rna_sobjs = epi_seurat_objs, promoter_length = c(5000, 1000), 
                          save_p2g_mat = T, gene_size = 3)

pdf("peak2gene_gs.pdf", width = 10, height = 10)
for(p in ps) {
  grid::grid.newpage()
  grid::grid.draw(p)
}
dev.off()

ps <- plot_peak2gene_link(proj_epi, genes = markerGenes, outdir = "peak2gene", use_groups = tissues, groupby = "tissue", 
                          use_matrix = "GeneIntegrationMatrix", rna_sobjs = epi_seurat_objs, promoter_length = c(5000, 1000), 
                          save_p2g_mat = T, gene_size = 3)

pdf("peak2gene.pdf", width = 10, height = 10)
for(p in ps) {
  grid::grid.newpage()
  grid::grid.draw(p)
}
dev.off()

# for(tissue in tissues) {
#   proj_one_group_path <- paste("peak2gene", gsub(" ", "_", tissue), sep = "/")
#   
#   seRNA <- readRDS(paste(proj_one_group_path, "Peak2GeneLinks/seRNA-Group-KNN.rds", sep = "/"))
#   seATAC <- readRDS(paste(proj_one_group_path, "Peak2GeneLinks/seATAC-Group-KNN.rds", sep = "/"))
#   
#   elementMetadata(seATAC)[["name"]] <- rowRanges(seATAC) %>% as.character
#   
#   colnames(seRNA) <- paste0("pseudo_cell_", seq_len( ncol(seRNA) ) )
#   colnames(seATAC) <- paste0("pseudo_cell_", seq_len( ncol(seATAC) ) )
#   
#   sRNA <- summarizedExperiment2Seurat(seRNA, assay = "RawRNA", do_norm = T, rename_assay = "RNA")
#   sATAC <- summarizedExperiment2Seurat(seATAC, assay = "RawATAC", do_norm = T, rename_assay = "ATAC")
#   
#   saveRDS(sRNA, file = paste(proj_one_group_path, "Peak2GeneLinks/seuratRNA-Group-KNN.rds", sep = "/") )
#   saveRDS(sATAC, file = paste(proj_one_group_path, "Peak2GeneLinks/seuratATAC-Group-KNN.rds", sep = "/") )
# }
# 
# for(tissue in tissues) {
#   proj_one_group_path <- paste("peak2gene", gsub(" ", "_", tissue), sep = "/")
#   
#   p2g <- readRDS(paste(proj_one_group_path, "peak2gene.rds", sep = "/"))
#   
#   p2g <- as.data.frame(p2g)
#   colnames(p2g)[c(1, 2)] <- c("peak", "gene")
#   
#   write.table(p2g, paste(proj_one_group_path, "peak2gene.txt", sep = "/"), sep='\t',
#               quote = F, row.names = F, col.names = T)
# }
# 
# for(tissue in tissues) {
#   cat(tissue, "\n")
#   proj_one_group_path <- paste("peak2gene", gsub(" ", "_", tissue), sep = "/")
#   
#   p2g_ht <- readRDS(paste(proj_one_group_path, "Peak2GeneLinks/p2g_heatmap_mat.rds", sep = "/"))
#   
#   p2g_filtered <- p2g_ht$Peak2GeneLinks
#   
#   write.table(p2g_filtered, file = paste(proj_one_group_path, "Peak2GeneLinks/filtered_peak2gene.txt", sep = "/"), sep = "\t", quote = F, col.names = NA)
# }


# TODO
# 1. pca irlba的错误是报在哪里的 √
# 2. 如何增加empty track √
# 3. 封装成函数与性能优化 √
# 4. 手动检查一下几个marker gene（特别是MUC13感觉有一点奇怪）√
# 5. 增加gene track √
# 6. 看一下cicero的画图实现（更加美观），viewpoint的实现 ×
# 7. check saveArchRProj的时间为什么会导致结果的区别 √
# 8. dropCells = T 的情况 √
# 9. 字体大小 √


