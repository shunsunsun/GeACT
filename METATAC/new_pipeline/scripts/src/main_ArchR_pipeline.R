#!/usr/bin/env r
suppressMessages({
  library(docopt)
})


# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [-s STAGE] [-t TISSUE] [--root ROOT] [--peaks PEAKS] [--merge_rna_stage RNASTAGE] [--cicero] [--fconns] [--threads THREADS]

-h --help           show help message
-s --stage STAGE    stage of sample [default: 11-14w]
-t --tissue TISSUE  tissue name [default: 01_stomach]
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
--peaks PEAKS       rds peaks file in results/peak_callong directory [default: normalPeaks.rds]
--merge_rna_stage RNASTAGE  RNA stage to use for integration and label transfer [default: 19-22w]
--cicero            run cicero pipeline [default: FALSE]
--fconns            force rerun conns step [default: FALSE]
--threads THREADS   Cicero threads [default: 16]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

if (opt$cicero %in% c("T", "F", "TRUE", "FALSE")){
  opt$cicero <- as.logical(opt$cicero)
} else{
  stop("--cicero not supported")
}

if (opt$fconns %in% c("T", "F", "TRUE", "FALSE")){
  opt$fconns <- as.logical(opt$fconns)
} else{
  stop("--fconns not supported")
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
  library(networkD3)
  library(plyr)
})

# envs -------------------------------------------------------------------------
set.seed(0)
stage <- opt$stage
tissue <- opt$tissue
runCicero <- opt$cicero
force_conns <- opt$fconns
root <- opt$root
peaks_file <- opt$peaks
rna_stage <- opt$merge_rna_stage
threads <- as.integer(opt$threads)

.data <- "data"
utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
expr_file <- paste(root, .data, stage, tissue, "RNA", rna_stage, "expr_RNA.rds", sep = "/")

organ_wd <- paste(root, .data, stage, tissue, sep = "/")
results_wd <- paste(organ_wd, "results", sep = "/")
ArchR_wd <- paste(results_wd, "ArchR", sep = '/')
Cicero_wd <- paste(results_wd, "Cicero", sep ="/")

source(utils_file)

# Open one PDF device to save all figures generated

setwd(results_wd)
pdf("pipeline.pdf", width = 12, height = 6, useDingbats = F, onefile = T)
pdf_dev <- dev.cur()

# 1 basic ArchR pipeline -------------------------------------------------------
setwd(ArchR_wd)
addArchRThreads(threads)
proj <- loadArchRProject(path = "ArchR_output", showLogo = F)

# load rna expr seurat object
cat(sprintf("Loading %s\n", expr_file))
expr <- readRDS(expr_file)
expr$tech <- "RNA"

# 1.1 based on tileMatrix ------------------------------------------------------
# dev.off()
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "tileIterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.8), 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  force = T
)

proj <- addTSNE(
  ArchRProj = proj, 
  reducedDims = "tileIterativeLSI", 
  name = "tileTSNE", 
  perplexity = 30, 
  force = T
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "tileIterativeLSI", 
  name = "tileUMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine", 
  force = T
 )

# add peak matrix --------------------------------------------------------------
peaks_file <- paste(root, .data, stage, tissue, "results/peak_calling", peaks_file, sep = "/")
# peaks_file <- paste(root, "meta/All_peaks_sorted_filt.narrowPeak", sep = "/")
peaks <- readRDS(peaks_file)
chrs <- paste0("chr", c(seq_len(22), "X"))
peaks <- GenomeInfoDb::keepSeqlevels(peaks, chrs, pruning.mode = "coarse")
proj <- addPeakSet(proj, peakSet = peaks, force = T)

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


# ArchR integration ------------------------------------------------------------
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "tileGeneIntegrationMatrix",
  reducedDims = "tileIterativeLSI",
  seRNA = expr,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "ident",
  nameCell = "predictedCell",
  nameGroup = "predictedIdent",
  nameScore = "predictedScore", 
  plotUMAP = T, 
  useImputation = F          # change this parameter to control whether use imputation during transfer
)

dev.set(pdf_dev)
# do the integration manually for visualization --------------------------------
# before integration
geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
geneScore <- summarizedExperiment2Seurat(se = geneScore, assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")
geneScore$tech <- "ATAC"

print(coembedding(geneScore, expr, features = VariableFeatures(expr), disc = "tech", assay1 = "ACTIVITY",
            subset = T, name1 = "ATAC activity", name2 = "RNA"))

# after integration
# data/label transfer
rD <- getReducedDims(ArchRProj = proj, reducedDims = "tileIterativeLSI", corCutOff = 0.75, dimsToUse = 1:30)
weight.reduction <- CreateDimReducObject(
  embeddings = rD, 
  key = "LSI_", 
  assay = NULL
)
geneScore_inte <- integrate(activity_atac = geneScore, expr_rna = expr, peak_matrix = peakMatrix, 
                            activity_slot = "ACTIVITY", weight.reduction = weight.reduction)

# visualization based on UMAP
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "predictedIdent"))

proj <- saveArchRProject(proj, load = T)

# Cicero pipeline --------------------------------------------------------------
if(runCicero){
  dir.create(Cicero_wd, recursive = T, showWarnings = F)
  setwd(Cicero_wd)
  
  cds <- peakMatrix2ciceroCDS(peakMatrix)
  input_cds <- cds$input_cds
  cicero_cds <- cds$cicero_cds
  
  # Run Cicero -------------------------------------------------------------------
  genome_size <- read.table(paste(root, "database", "txdb", "hg38.chrom22X.sizes", sep = "/"), sep='\t', header = F)
  
  if(file.exists("conns.txt.gz") && !force_conns){
    conns <- read.table(gzfile("../Cicero/conns.txt.gz"), sep = "\t", header = T)
  } else {
    conns <- run_cicero(cicero_cds, genome_size, sample_num = 100)
    write.table(conns, gzfile("conns.txt.gz"), sep = '\t', row.names = F, col.names = T, quote = F)
  }
  
  # Finding cis-Co-accessibility Networks (CCANS)
  CCAN_assigns <- generate_ccans(conns)
  write.table(CCAN_assigns, gzfile("CCAN_assigns.txt.gz"), sep='\t', row.names = T, col.names = T, quote = F)
  
  # Cicero gene activity scores
  # Add a column for the pData table indicating the gene if a peak is a promoter
  # Create a gene annotation set that only marks the transcription start sites of 
  # the genes. We use this as a proxy for promoters.
  gene_annotation_sub <- read.table(paste(root, "database", "human", "gene_position.hg38.v26.bed", sep = "/"),
                                    sep='\t', header = F, stringsAsFactors = F)
  chr_ids <- paste0("chr", c(1:22, "X"))
  colnames(gene_annotation_sub) <- c("chromosome", "start", "end", "gene", "ensemblid", "strand")
  gene_annotation_sub <- subset(gene_annotation_sub, is.element(gene_annotation_sub$chromosome, chr_ids))
  gene_annotation_sub[gene_annotation_sub$strand == "+", "start"] <- gene_annotation_sub[gene_annotation_sub$strand == "+", "start"] - 2000
  gene_annotation_sub[gene_annotation_sub$strand == "-", "end"] <- gene_annotation_sub[gene_annotation_sub$strand == "-", "end"] + 2000
  # gene_annotation_sub[gene_annotation_sub$strand == "+", "end"] <- gene_annotation_sub[gene_annotation_sub$strand == "+", "start"] + 1
  # gene_annotation_sub[gene_annotation_sub$strand == "-", "start"] <- gene_annotation_sub[gene_annotation_sub$strand == "-", "end"] - 1
  
  gene_id2name <- read.table(paste(root, "database", "GENCODE", "v26", "gene_ID2Name_fixed.txt", sep = "/"), header = F, stringsAsFactors = F)
  gene_id2name <- column_to_rownames(gene_id2name, "V1")
  colnames(gene_id2name) <- "gene"
  gene_annotation_sub$gene <- gene_id2name[gene_annotation_sub$ensemblid, ]
  
  input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
  
  write.table(fData(input_cds), "annotate.txt", sep='\t', row.names = F, col.names = T, quote = F)
  
  
  # Generate gene activity scores
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # remove any rows/columns with all zeroes
  unnorm_ga <- unnorm_ga[Matrix::rowSums(unnorm_ga) != 0, 
                         Matrix::colSums(unnorm_ga) != 0]
  write.table(as.matrix(unnorm_ga), gzfile("cicero_gene_activities_unnorm.tsv.gz"),
              sep = '\t',row.names = T,col.names = T, quote = F)
  
  # normalize
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  norm.ga <- normalize_gene_activities(unnorm_ga, num_genes)
  write.table(as.matrix(norm.ga),gzfile("cicero_gene_activities.tsv.gz"),
              sep = '\t',row.names = T,col.names = T, quote = F)
  
  # 2.2 Coembedding of cicero gene activity score
  activity.matrix <- read.table(gzfile("cicero_gene_activities_unnorm.tsv.gz"),
                                sep='\t', header = T, stringsAsFactors = F, check.names = F, comment.char = "")
  activity_atac <- CreateSeuratObject(counts = as.sparse(activity.matrix), assay = "ACTIVITY", 
                                      project = "Geact_ATAC")
  activity_atac$tech <- "ATAC"
  print(coembedding(activity_atac, expr, features = VariableFeatures(expr), disc = "tech", assay1 = "ACTIVITY", subset = T, 
              name1 = "Cicero ATAC activity", name2 = "RNA"))
  cicero_inte <- integrate(activity_atac = activity_atac, expr_rna = expr, peak_matrix = peakMatrix, 
                           activity_slot = "ACTIVITY")
  
  # compare between geneScore and cicero results
  lev <- unique(expr[[]][["ident"]])
  p <- cmp_transfer(geneScore_inte$celltype.predictions$predicted.id, cicero_inte$celltype.predictions$predicted.id, lev = lev)
  saveNetwork(p, "geneScore_cicero.html")
  
  # save cicero prediction
  proj$ciceroIdent <- cicero_inte$celltype.predictions["predicted.id"][getCellNames(proj), ]
}


# Generate Output Info ---------------------------------------------------------
setwd(results_wd)
# peak and tile matrix
dir.create("peak_matrix", showWarnings = F)
writeMM(assay(peakMatrix), file = "peak_matrix/peak_matrix.mtx")
write.table(colnames(peakMatrix), file = "peak_matrix/cellname.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(rownames(peakMatrix), file = "peak_matrix/peakinfo.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# meta table
cellMeta <- getCellColData(proj)[, c("seqID", "tissue", "samplingPos", "plate", "individual", "PassQC", "predictedIdent", "Reads", "Aligned_ratio",
                                "nFrags", "FRIP", "mito_ratio", "BlacklistRatio", "DoubletScore")]
colnames(cellMeta) <- c("seqID", "tissue", "samplingPos", "plate", "individual", "QC", "ident", "cleanReads", "mpRatio",
                        "nFrags", "FRIP", "MitoRatio", "BlacklistRatio", "DoubletScore")
cellMeta$nPeak <- colSums(assay(peakMatrix) > 0)[rownames(cellMeta)]
cellMeta$stage <- stage
cellMeta$species <- "human"
cellMeta$QC <- as.logical(cellMeta$QC)
cellMeta$group <- ident2clgrp(cellMeta$ident)
proj$group <- cellMeta[, "group", drop = F][getCellNames(proj), ]
proj$nPeak <- cellMeta[, "nPeak", drop = F][getCellNames(proj), ]

# add umap/tsne embedding
peakTSNE <- getEmbedding(proj, embedding = "peakTSNE")
peakUMAP <- getEmbedding(proj, embedding = "peakUMAP")

colnames(peakTSNE) <- c("tSNE_1", "tSNE_2")
colnames(peakUMAP) <- c("UMAP_1", "UMAP_2")

cellMeta <- cbind(cellMeta, peakTSNE)
cellMeta <- cbind(cellMeta, peakUMAP)

# add Cicero results
if(runCicero){
  cellMeta$ciceroIdent <- cicero_inte$celltype.predictions["predicted.id"][rownames(cellMeta), ]
}

write.table(cellMeta, file = "filtered_cellMeta_internal.txt", sep = "\t", quote = F, col.names = NA)

rownames(cellMeta) <- gsub("^.*#", "", rownames(cellMeta))
rownames(geneScore_inte$tsne.data)[geneScore_inte$tsne.data$tech == "atac"] <- rownames(geneScore_inte$tsne.data)[geneScore_inte$tsne.data$tech == "atac"] %>% gsub("^.*#", "", .)


if (stage == "19-22w") {
  rownames(cellMeta) <- gsub("^(.*?_)", "\\1A_", rownames(cellMeta)) 
  rownames(geneScore_inte$tsne.data)[geneScore_inte$tsne.data$tech == "atac"] <- rownames(geneScore_inte$tsne.data)[geneScore_inte$tsne.data$tech == "atac"] %>% gsub("^(.*?_)", "\\1A_", .)
} else{
  rownames(cellMeta) <- gsub("^(.*?_)", "\\1B_", rownames(cellMeta))
  rownames(geneScore_inte$tsne.data)[geneScore_inte$tsne.data$tech == "atac"] <- rownames(geneScore_inte$tsne.data)[geneScore_inte$tsne.data$tech == "atac"] %>% gsub("^(.*?_)", "\\1B_", .)
}

write.table(cellMeta, file = "filtered_cellMeta.txt", sep = "\t", quote = F, col.names = NA)

# write integration dimension reduction results 
write.table(geneScore_inte$tsne.data, file = "integration_dimreduc.txt", sep = "\t", quote = F, col.names = NA)

# close pdf device
dev.off()

# save final ArchR project
saveArchRProject(proj, load = F)

# dimension reduction visualization --------------------------------------------
pdf("dimRed.pdf", width = 10, height = 10)
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "seqID", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "samplingPos", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "plate", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "individual", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "predictedIdent", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
if (runCicero)
  print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "ciceroIdent", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "nFrags", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "FRIP", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "mito_ratio", plotAs = "points", size = 1.5))
print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "DoubletScore", plotAs = "points", size = 1.5))
dev.off()

# marker gene visualization ----------------------------------------------------
# generate geneScore again to get group info
geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix") %>% 
  summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")

dir.create("markers", showWarnings = F)

pdf("markers/simple_markers1.pdf", width = 10, height = 10)
simpleMarkersPlot1(proj, geneScore)
dev.off()

pdf("markers/simple_markers2.pdf", width = 10, height = 10)
simpleMarkersPlot2(proj, geneScore, tissue = tissue, root = root)
dev.off()

pdf("markers/signature_genes.pdf", width = 10, height = 10)
signatureGenesPlot(proj, geneScore, tissue = tissue, root = root)
dev.off()
