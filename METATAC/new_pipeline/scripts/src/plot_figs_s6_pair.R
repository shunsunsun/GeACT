#!/usr/bin/env r
suppressPackageStartupMessages({
  library(docopt)
})

# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [--root ROOT] [--tissue TISSUE]

-h --help           show help message
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
-t --tissue TISSUE  tissue to analysis [default: 01_stomach]
"

opt <- docopt(doc)

root <- opt$root
tissue <- opt$tissue

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

setwd(paste("/data/Lab/otherwork/GeACT/ATAC/data/paired", tissue, sep = "/"))
set.seed(seed = 0)

stages <- c("11-14w", "19-22w")
rna_files <- paste(root, "data", stages, tissue, "RNA", stages, "expr_RNA.rds", sep = "/")
if(!file.exists(rna_files[1])) {
  rna_files[1] <- paste(root, "data", "11-14w", tissue, "RNA", "19-22w", "expr_RNA.rds", sep = "/")
}
names(rna_files) <- stages
names(stages) <- stages

seurat_objs <- lapply(stages, function(stage) {
  cat(stage, "\n")
  sobj <- readRDS(rna_files[stage])
  sobj$group <- ident2clgrp(sobj$ident)
  sobj$stage <- stage
  sobj
})

proj <- loadArchRProject("ArchR_output/", showLogo = F)

ps <- plot_peak2gene_link(proj, genes = "CLDN18", outdir = "peak2gene", groupby = "stage", 
                          use_matrix = "GeneIntegrationMatrix", promoter_length = c(5000, 10000), 
                          save_p2g_mat = T, rna_sobjs = seurat_objs)

# grid::grid.draw(ps[[1]])
