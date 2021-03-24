#!/usr/bin/env r
suppressMessages({
  library(docopt)
})

# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [-t TISSUE] [--root ROOT] [-u UP] [--threads THREADS]

-h --help           show help message
-t --tissue TISSUE  tissue name [default: 01_stomach]
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
-u --unionPeaks UP  union peaks name to use [default: unionPeaks.rds]
--threads THREADS   Cicero threads [default: 16]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

suppressMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
})

# envs -------------------------------------------------------------------------
set.seed(0)
tissue <- opt$tissue
root <- opt$root
threads <- as.integer(opt$threads)

results_wd <- paste(root, "data", "paired", tissue, sep = "/")

stages <- c("19-22w", "11-14w")

frag_files <- sapply(stages, function(stage){paste(root, "data", stage, tissue, "frag_and_meta", "merge_human_frag_decon.bed.gz", sep = "/")})

unionPeaks_file <- paste(root, "data", "paired", tissue, opt$unionPeaks, sep = "/")

# load data
unionPeaks <- readRDS(unionPeaks_file)

cell_metas <- lapply(stages, function(stage){
  cell_meta_file <- paste(root, "data", stage, tissue, "results", "filtered_cellMeta_internal.txt", sep = "/")
  cell_meta <- read.table(cell_meta_file, header = T, sep = "\t", quote = "", comment.char = "") %>% column_to_rownames(var = "X")
})
cell_meta <- rbind(cell_metas[[1]], cell_metas[[2]])

geneAnnotation <- readRDS(paste(root, "database", "annotation", "geneAnnotation.rds", sep = "/"))
genomeAnnotation <- readRDS(paste(root, "database", "annotation", "genomeAnnotation.rds", sep = "/"))

# ArchR pipeline
setwd(results_wd)
addArchRThreads(threads)

ArrowFiles <- createArrowFiles(
  inputFiles = frag_files,
  sampleNames = c(paste(stages, tissue, sep = '_')),
  minTSS = 0,
  minFrags = 0,
  maxFrags = Inf,
  addTileMat = F,
  addGeneScoreMat = F,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  force = T,
  subThreading = F
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR_output",
  copyArrows = F,
  geneAnnotation = geneAnnotation, 
  genomeAnnotation = genomeAnnotation, 
  showLogo = F
)

# filtering
cell_names <- rownames(cell_meta)
proj <- proj[cell_names]

for (i in seq_along(cell_meta)) {
  proj <- addCellColData(proj, data = cell_meta[[i]], name = names(cell_meta[i]), cells = rownames(cell_meta), force = T)
}

# add peakset and peak matrix
chrs <- paste0("chr", c(seq_len(22), "X"))
unionPeaks <- GenomeInfoDb::keepSeqlevels(unionPeaks, chrs, pruning.mode = "coarse")
proj <- addPeakSet(proj, peakSet = unionPeaks, force = T)

# counting
proj <- addPeakMatrix(proj, ceiling = 100, force = T)
cat(getAvailableMatrices(proj), "\n")

peakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
peaks <- rowRanges(peakMatrix)
rownames(peakMatrix) <- paste(as.character(seqnames(peaks)), start(peaks), end(peaks), sep = "_")

# add embedding based on peaks
proj <- addIterativeLSI(proj, useMatrix = "PeakMatrix", 
                        iterations = 1, name = "peakLSI", varFeatures = 50000, force = T)

proj <- addTSNE(
  ArchRProj = proj, 
  reducedDims = "peakLSI", 
  name = "peakTSNE", 
  perplexity = 30, 
  force = T
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "peakLSI", 
  name = "peakUMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine", 
  force = T
)

# add paired tSNE to cell meta
peakTSNE <- getEmbedding(proj, embedding = "peakTSNE")
peakUMAP <- getEmbedding(proj, embedding = "peakUMAP")

colnames(peakTSNE) <- c("ali_tSNE_1", "ali_tSNE_2")
colnames(peakUMAP) <- c("ali_UMAP_1", "ali_UMAP_2")

cell_meta <- cbind(cell_meta, peakTSNE)
cell_meta <- cbind(cell_meta, peakUMAP)

write.table(cell_meta, file = "cellMeta_internal.txt", sep = "\t", quote = F, col.names = NA)

rownames(cell_meta)[cell_meta$stage == "19-22w"] <- gsub("#(.*?_)", "#\\1A_", rownames(cell_meta)[cell_meta$stage == "19-22w"]) 
rownames(cell_meta)[cell_meta$stage == "11-14w"] <- gsub("#(.*?_)", "#\\1B_", rownames(cell_meta)[cell_meta$stage == "11-14w"])

rownames(cell_meta) <- gsub("^.*#", "", rownames(cell_meta))

write.table(cell_meta, file = "cellMeta.txt", sep = "\t", quote = F, col.names = NA)

# save embeddings
pdf("paired_stage.pdf", width = 10, height = 10)
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "stage", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "ident", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", size = 1.5))
dev.off()
