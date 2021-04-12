#!/usr/bin/env r
suppressMessages({
  library(docopt)
})

# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [-t TISSUE] [--root ROOT] [-u UP] [--threads THREADS]

-h --help           show help message
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
-u --unionPeaks UP  union peaks name to use in data/all [default: unionPeaks.rds]
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

root <- opt$root
threads <- as.integer(opt$threads)
unionPeaks_file <- opt$unionPeaks

# env
results_wd <- paste(root, "data", "all", sep = "/")



setwd(results_wd)
addArchRThreads(threads)

geneAnnotation <- readRDS(paste(root, "database", "annotation", "geneAnnotation.rds", sep = "/"))
genomeAnnotation <- readRDS(paste(root, "database", "annotation", "genomeAnnotation.rds", sep = "/"))

# config input files
stages <- c("19-22w", "11-14w")
stage_tissues <- list("19-22w"=c("01_stomach", "02_small_intestine", "03_kidney", "04_lung", "05_pancreas", "06_spleen", 
                            "07_testis", "08_bladder", "09_bone_marrow", "11_diaphragm", "12_esophagus", "14_large_intestine", "15_liver", "16_ovary", "17_thymus"),
                     "11-14w"=c("01_stomach", "02_small_intestine", "03_kidney", "04_lung", "05_pancreas", "08_bladder", "09_bone_marrow", "10_bronchus", "12_esophagus", "15_liver", "16_ovary", "17_thymus"))
pair_tissues <- intersect(stage_tissues[[1]], stage_tissues[[2]])

frag_files <- c()
cell_meta_files <- c()
sample_name <- c()

for(stage in stages){
  for (tissue in stage_tissues[[stage]]){
    sample_name <- c(sample_name, paste(stage, tissue, sep = "_"))
    frag_files <- c(frag_files, paste(root, "data", stage, tissue, "frag_and_meta/merge_human_frag_decon.bed.gz", sep = "/"))
    cell_meta_files <- c(cell_meta_files, paste(root, "data", stage, tissue, "results/filtered_cellMeta_internal.txt", sep = "/"))
  }
}

# load cell meta
cell_meta <- NULL
for(file_name in cell_meta_files){
  tmp_meta <- read_tsv(file_name, col_names = T, col_types = cols(group = col_character()), quote = "") %>% column_to_rownames(var = "X1")
  cell_meta <- rbind(cell_meta, tmp_meta)
}

# load paired meta
pair_cell_meta <- NULL
for(tissue in pair_tissues){
  file_name <- paste(root, "data/paired", tissue, "cellMeta_internal.txt", sep = "/")
  tmp_meta <- read_tsv(file_name, col_names = T, col_types = cols(group = col_character()), quote = "") %>% column_to_rownames(var = "X1")
  pair_cell_meta <- rbind(pair_cell_meta, tmp_meta)
}

cell_meta <- merge(cell_meta, pair_cell_meta[c("pair_tSNE_1", "pair_tSNE_2", "pair_UMAP_1", "pair_UMAP_2")], by = 0, all.x = T) %>% column_to_rownames("Row.names")

# load peaks
unionPeaks <- readRDS(unionPeaks_file)


# ArchR pipeline
ArrowFiles <- createArrowFiles(
  inputFiles = frag_files,
  sampleNames = sample_name,
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

saveArchRProject(proj, load = F)
# proj <- loadArchRProject("ArchR_output")

# add paired tSNE to cell meta
peakTSNE <- getEmbedding(proj, embedding = "peakTSNE")
peakUMAP <- getEmbedding(proj, embedding = "peakUMAP")

colnames(peakTSNE) <- c("all_tSNE_1", "all_tSNE_2")
colnames(peakUMAP) <- c("all_UMAP_1", "all_UMAP_2")

cell_meta <- cbind(cell_meta, peakTSNE)
cell_meta <- cbind(cell_meta, peakUMAP)

write.table(cell_meta, file = "cellMeta_internal.txt", sep = "\t", quote = F, col.names = NA)

rownames(cell_meta)[cell_meta$stage == "19-22w"] <- gsub("#(.*?_)", "#\\1A_", rownames(cell_meta)[cell_meta$stage == "19-22w"]) 
rownames(cell_meta)[cell_meta$stage == "11-14w"] <- gsub("#(.*?_)", "#\\1B_", rownames(cell_meta)[cell_meta$stage == "11-14w"])

rownames(cell_meta) <- gsub("^.*#", "", rownames(cell_meta))

write.table(cell_meta, file = "cellMeta.txt", sep = "\t", quote = F, col.names = NA)

# save embeddings
pdf("merge_all.pdf", width = 10, height = 20)
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "stage", size = 1))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "ident", size = 1))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", size = 1))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tissue", size = 1))
dev.off()
