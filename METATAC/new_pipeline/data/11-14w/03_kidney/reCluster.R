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
  library(ArchR)
})

root <- opt$root
stage <- "11-14w"

setwd(paste(root, "data/11-14w/03_kidney/results", sep = "/"))

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

proj <- loadArchRProject("ArchR/ArchR_output/", showLogo = F)

cellMeta <- getCellColData(proj)[, c("seqID", "tissue", "samplingPos", "plate", "individual", "PassQC", "predictedIdent", "Reads", "Aligned_ratio",
                                     "nFrags", "FRIP", "mito_ratio", "BlacklistRatio", "DoubletScore", "nPeak", "group")]
colnames(cellMeta) <- c("seqID", "tissue", "samplingPos", "plate", "individual", "QC", "ident", "cleanReads", "mpRatio",
                        "nFrags", "FRIP", "MitoRatio", "BlacklistRatio", "DoubletScore", "nPeak", "group")

cellMeta$tuned_group <- cellMeta$group

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
