library(ArchR)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(gridExtra)
library(Seurat)
library(cicero)
# library(reticulate)
# use_condaenv("macs2")

# envs -------------------------------------------------------------------------
age <- '19-22w'
tissue <- '02_small_intestine'
root <- "/data/Lab/otherwork/GeACT/ATAC"
.data <- "data"

organ_wd <- paste(root, .data, age, tissue, sep = '/')
ArchR_wd <- paste(root, .data, age, tissue, "results/ArchR", sep = '/')
Cicero_path <- paste(root, .data, age, tissue, "results/Cicero", sep = '/')
setwd(ArchR_wd)


# utils ------------------------------------------------------------------------
my_clip <- function(x, lb, ub){
  pmax(lb, pmin(x, ub))
}

change_assay <- function(sobj, assayin = "RNA", assayout = "ACTIVITY"){
  sobj[[assayout]] <- sobj[[assayin]]
  DefaultAssay(sobj) <- assayout
  sobj[[assayin]] <- NULL
  
  return(sobj)
}

coembedding <- function(sobj1, sobj2, features = NULL, disc = "tech", renorm = F,
                        assay1 = "RNA", assay2 = "RNA", subset = F){
  if(assay1 != "RNA") sobj1 <- change_assay(sobj1, assay1, "RNA")
  if(assay2 != "RNA") sobj2 <- change_assay(sobj2, assay2, "RNA")
  if (subset) {
    genes_both <- intersect(rownames(sobj1), rownames(sobj2))
    sobj1 <- sobj1[genes_both, ]
    sobj2 <- sobj2[genes_both, ]
  }
  merge_obj <- merge(sobj1, sobj2)
  cat(sprintf("Merge obj of (%d, %d) and (%d, %d) into (%d, %d)\n", 
                dim(sobj1)[1], dim(sobj1)[2],
                dim(sobj2)[1], dim(sobj2)[2],
                dim(merge_obj)[1], dim(merge_obj)[2]))

  if (renorm) merge_obj <- NormalizeData(merge_obj, normalization.method = "LogNormalize")
  features <- if (is.null(features)) VariableFeatures(merge_obj) else features
  merge_obj <- ScaleData(merge_obj, features = features)
  merge_obj <- RunPCA(merge_obj, features = features, npcs = 50, verbose = F)
  merge_obj <- RunTSNE(merge_obj, dims = 1:30)
  DimPlot(merge_obj, group.by = disc)
}

summarizedExperiment2Seurat <- function(se, assay = "RNA", project = "Project", do_norm = F){
  gene_meta <- as.data.frame(elementMetadata(se)) %>% column_to_rownames("name")
  gene_meta <- gene_meta[colnames(gene_meta) != "idx"]
  rownames(se) <- rownames(gene_meta)
  
  cell_meta <- as.data.frame(colData(se))
  
  se_assay <- assay(se)
  
  sobj <- CreateSeuratObject(counts = se_assay, project = project, meta.data = cell_meta, assay = assay)
  sobj[[assay]] <- AddMetaData(sobj[[assay]], gene_meta)
  
  if(do_norm) sobj <- NormalizeData(sobj, normalization.method = "LogNormalize")
  return(sobj)
}

integrate <- function(activity_atac, expr_rna, peak_matrix, upper_bound = 5, activity_slot = "ACTIVITY", 
                      weight.reduction = NULL){
  # atac object setup
  activity_atac[["ATAC"]] <- CreateAssayObject(counts = assay(peak_matrix))
  activity_atac$tech <- "atac"
  
  DefaultAssay(activity_atac) <- activity_slot
  activity_atac <- NormalizeData(activity_atac)
  # print(ggplot(data = rowMeans(activity_atac, slot = "counts") %>% my_clip(lb = 0.1, ub = upper_bound) %>% as_tibble, aes(x=value)) + geom_histogram(bins = 20))
  # activity_atac <- FindVariableFeatures(activity_atac, num.bin = 20, mean.cutoff = c(0.1, upper_bound))
  activity_atac <- ScaleData(activity_atac, features = rownames(activity_atac))
  
  DefaultAssay(activity_atac) <- "ATAC"
  VariableFeatures(activity_atac) <- names(which(Matrix::rowSums(activity_atac) > ncol(activity_atac)))
  activity_atac <- RunLSI(activity_atac, n = 50, scale.max = NULL)
  activity_atac <- RunUMAP(activity_atac, reduction = "lsi", dims = 1:50)
  # expr_atac <- NormalizeData(expr_atac)
  # expr_atac <- ScaleData(expr_atac, features = rownames(expr_atac[['ATAC']]))
  
  # rne object setup
  expr_rna <- NormalizeData(expr_rna)
  expr_rna[['ident']] <- Idents(expr_rna)
  expr_rna$tech <- "rna"
  genes.use <- intersect(rownames(expr_rna), rownames(activity_atac[[activity_slot]]))
  expr_rna_sub <- expr_rna[genes.use, ]
  expr_rna_sub <- FindVariableFeatures(expr_rna_sub, nfeatures = 2000)
  
  # find anchors and transfer data
  transfer.anchors <- FindTransferAnchors(
    reference = expr_rna_sub, query = activity_atac,
    features = VariableFeatures(object = expr_rna_sub),  #TOTEST 测试是否需要intersect
    normalization.method = "LogNormalize",
    reference.assay = "RNA", query.assay = activity_slot, reduction = "cca")
  
  if (is.null(weight.reduction)) weight.reduction <- activity_atac[["lsi"]]
  
  celltype.predictions <- TransferData(
    anchorset = transfer.anchors, refdata = expr_rna_sub$ident,
    weight.reduction = weight.reduction)
  
  activity_atac <- AddMetaData(activity_atac, metadata = celltype.predictions$predicted.id, "predicted.id")
  activity_atac <- AddMetaData(activity_atac, metadata = celltype.predictions$prediction.score.max, "prediction.score.max")
  
  # to make the colors match
  activity_atac$predicted.id <- factor(activity_atac$predicted.id, levels = levels(expr_rna_sub))
  Idents(activity_atac) <- celltype.predictions$predicted.id
  # write.table(celltype.predictions,"../Seurat_integration/celltype_predictions.txt",
  #             sep='\t',quote = F,row.names = T,col.names = T)
  
  # Co-embedding
  refdata <- GetAssayData(expr_rna_sub, assay = "RNA", slot = "data")
  
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                             weight.reduction = weight.reduction)
  cat(paste0("Tranfer dim: ", dim(imputation)[1], " ", dim(imputation)[2]))
  
  activity_atac[["RNA"]] <- imputation
  # saveRDS(expr.atac, file = "../Seurat_integration/Seurat_expr_atac_integration.rds")
  coembed <- merge(x = expr_rna_sub, y = activity_atac)
  
  # Finally, we run PCA and TSNE on this combined object, to visualize the co-embedding of both
  # datasets
  coembed <- ScaleData(coembed, features = genes.use, do.scale = F, verbose = F)
  coembed <- RunPCA(coembed, features = genes.use, verbose = F)
  coembed <- RunTSNE(coembed, dims = 1:30)
  coembed$celltype <- ifelse(!is.na(coembed$ident), coembed$ident, coembed$predicted.id)
  # saveRDS(coembed, file = "../Seurat_integration/Seurat_expr_coembed.rds")
  
  p1 <- DimPlot(coembed, group.by = "tech")
  p2 <- DimPlot(coembed, group.by = "ident", label = T, repel = T)
  print(p1 + p2)
  
  list(imputation = imputation, celltype.predictions = celltype.predictions)
}


# set up ArchR------------------------------------------------------------------
addArchRThreads(threads = 32)
# addArchRGenome("hg38")

proj <- loadArchRProject("ArchR_output", showLogo = F)
expr <- readRDS("../../RNA/expr_RNA.rds")
expr$tech <- "RNA"

# 1 New pipeline ---------------------------------------------------------------
# generate peak matrix
peaks_file <- paste(root, .data, age, tissue, "results/peak_calling/macs2/organ_peaks_filtered.narrowPeak", sep = "/")
# peaks_file <- paste(root, "meta/All_peaks_sorted_filt.narrowPeak", sep = "/")
peaks <- rtracklayer::import(peaks_file)
chrs <- paste0("chr", c(seq_len(22), "X"))
peaks <- GenomeInfoDb::keepSeqlevels(peaks, chrs, pruning.mode = "coarse")
proj <- addPeakSet(proj, peakSet = peaks, force = T)

# counting
proj <- addPeakMatrix(proj, ceiling = 100)
getAvailableMatrices(proj)

# 1.0 check genename between ATAC gene score and RNA ---------------------------
# these steps make nomore sense, due to seurat substitution of "_" with "-"
# genes_in_ATAC <- (getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = F) %>% elementMetadata)$name
# genes_in_RNA <- rownames(expr)
# 
# genes_both <- intersect(genes_in_ATAC, genes_in_RNA)
# 
# # subset of genes in expr (note normalization have been done for expr object, please don't do normalization in downstream analysis)
# expr_sub <- expr[genes_both, ]
# expr_sub <- FindVariableFeatures(expr_sub, nfeatures = 2000)

# 1.1 GeneScore matrix ---------------------------------------------------------
# 1.1.1 directly visualize -----------------------------------------------------
# without subsetting to shared gene list
geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
geneScore <- summarizedExperiment2Seurat(se = geneScore, assay = "ACTIVITY")
geneScore$tech <- "ATAC"
coembedding(geneScore, expr, features = VariableFeatures(expr), disc = "tech", assay1 = "ACTIVITY")
# after subsetting
coembedding(geneScore, expr, features = VariableFeatures(expr), disc = "tech", assay1 = "ACTIVITY",
            subset = T)

# 1.1.2 Seurat integration of genescore matrix(different levels) ---------------
# load peakMatrix(need LSI reduction for integration visualization)
# visual check (head and LSI dimension reduction and visualization)
peakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rownames(peakMatrix) <- rowRanges(peakMatrix)$name

dev.off()
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
  perplexity = 30
)

proj <- addIterativeLSI(proj, useMatrix = "PeakMatrix", iterations = 1, name = "peakLSI")

proj <- addTSNE(
  ArchRProj = proj, 
  reducedDims = "peakLSI", 
  name = "peakTSNE", 
  perplexity = 30
)

plotEmbedding(proj, embedding = "tileTSNE", imputeWeights = NULL)
plotEmbedding(proj, embedding = "peakTSNE", imputeWeights = NULL)

# data/label transfer
rD <- getReducedDims(ArchRProj = proj, reducedDims = "tileIterativeLSI", corCutOff = 0.75, dimsToUse = 1:30)
weight.reduction <- CreateDimReducObject(
  embeddings = rD, 
  key = "LSI_", 
  assay = NULL
)
geneScore_inte <- integrate(activity_atac = geneScore, expr_rna = expr, peak_matrix = peakMatrix, 
                            upper_bound = 3, activity_slot = "ACTIVITY", weight.reduction = weight.reduction)
# geneScore_inte <- integrate(activity_atac = geneScore, expr_rna = expr, peak_matrix = peakMatrix, 
#                             upper_bound = 3, activity_slot = "ACTIVITY")


# 1.2 GeneIntegrationMatrix ----------------------------------------------------
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
#TODO 为什么UMAP输出少了一张图？
#TODO warning message?

# 1.2.1 jointCCA visualization(inside ArchR) -----------------------------------
# load CCA results directly, and plot tsne based upon cca:1-30
jointCCA <- readRDS("./ArchR_output/RNAIntegration/tileGeneIntegrationMatrix/Save-Block1-JointCCA.rds")
set.seed(1)
UMAPParams <- list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE)
UMAPParams$X <- as.data.frame(jointCCA[, grep("CC_", colnames(jointCCA))])
UMAPParams$ret_nn <- FALSE
UMAPParams$ret_model <- FALSE
UMAPParams$n_threads <- 1
uwotUmap <- do.call(uwot::umap, UMAPParams)
jointCCA$UMAP1 <- uwotUmap[, 1]
jointCCA$UMAP2 <- uwotUmap[, 2]
ggplot(data = as.data.frame(jointCCA), aes(x = UMAP1, y = UMAP2, colour = Assay)) + geom_point(size = 0.3) + theme_ArchR()


# 1.2.2 coembedding of GeneIntegrationMatrix -----------------------------------
geneIntegrationMatrix <- getMatrixFromProject(proj, useMatrix = "tileGeneIntegrationMatrix", verbose = F)
geneIntegrationMatrix <- summarizedExperiment2Seurat(geneIntegrationMatrix, assay = "ACTIVITY")
geneIntegrationMatrix$tech <- "ATAC"

# no subsetting
coembedding(geneIntegrationMatrix, expr, features = VariableFeatures(expr), disc = "tech", assay1 = "ACTIVITY")

# subset
geneIntegrationMatrix_sub <- geneIntegrationMatrix[genes_both, ]
coembedding(geneIntegrationMatrix_sub, expr_sub, features = VariableFeatures(expr_sub), disc = "tech", assay1 = "ACTIVITY")

# 1.3 compare self implemented and archr labels --------------------------------
# Note: if implemented with DimReduction object from ArchR, it obtains exactly the same results compared to ArchR
archr_result <- getCellColData(proj)[, c("predictedIdent", "predictedScore")]
self_result <- geneScore_inte$celltype.predictions[, "predicted.id", drop = F]
cmp <- merge(archr_result, self_result, by.x = 0, by.y = 0)
sum(cmp$predictedIdent == cmp$predicted.id)



# 2 Old pipeline ---------------------------------------------------------------
# 2.1 Cicero pipeline ----------------------------------------------------------
force_cicero <- F
dir.create(Cicero_path, recursive = T, showWarnings = F)
setwd(Cicero_path)

peakMatrix2ciceroCDS <- function(peakMatrix){
  indata <- assay(peakMatrix)
  # binarize the matrix
  indata@x[indata@x > 0] <- 1
  
  # format cell info
  cellinfo <- as.data.frame(colData(peakMatrix))
  
  # format peak info
  peaks <- rowRanges(peakMatrix)
  peakinfo <- data.frame(row.names = paste(as.character(seqnames(peaks)), start(peaks), end(peaks), sep = "_"), 
                         name = peaks$name)
  
  rownames(peakMatrix) <- rownames(peakinfo)
  
  rownames(indata) <- rownames(peakinfo)
  colnames(indata) <- rownames(cellinfo)
  
  # make CDS
  input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                  cell_metadata = cellinfo,
                                                  gene_metadata = peakinfo))
  
  #Ensure there are no peaks included with zero reads
  dim(input_cds)
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]
  dim(input_cds)
  
  #### dimensionality reduction by UMAP
  set.seed(2020)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
  plot_cells(input_cds)
  
  #### access the UMAP coordinates
  umap_coords <- reducedDims(input_cds)$UMAP
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
}

cicero_cds <- peakMatrix2ciceroCDS(peakMatrix)

# Run Cicero -------------------------------------------------------------------
genome_size <- read.table(paste(root, "database", "txdb", "hg38.chrom22X.sizes", sep = "/"), sep='\t', header = F)

if(file.exists("conns.txt.gz") && !force_cicero){
  conns <- read.table(gzfile("../Cicero/conns.txt.gz"), sep = "\t", header = T)
} else {
  conns <- run_cicero(cicero_cds, genome_size, sample_num = 100)
  write.table(conns, gzfile("../Cicero/conns.txt.gz"), sep = '\t', row.names = F, col.names = T, quote = F)
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
coembedding(activity_atac, expr, features = VariableFeatures(expr), disc = "tech", assay1 = "ACTIVITY", subset = T)
cicero_inte <- integrate(activity_atac = activity_atac, expr_rna = expr, peak_matrix = peakMatrix, upper_bound = 5,
          activity_slot = "ACTIVITY")


# Compare ArchR and Cicero results
archr_result <- getCellColData(proj)[, c("predictedIdent", "predictedScore")]
self_result <- cicero_inte$celltype.predictions[, "predicted.id", drop = F]
cmp <- merge(archr_result, self_result, by.x = 0, by.y = 0)
sum(cmp$predictedIdent == cmp$predicted.id)


# 3 Generation of downstream analysis data


# 4 Other experiments



# utils and experiments --------------------------------------------------------

my_addIterativeLSI <- function (ArchRProj = NULL, useMatrix = "TileMatrix", name = "IterativeLSI", 
          iterations = 2, clusterParams = list(resolution = c(2), 
                                               sampleCells = 10000, maxClusters = 6, n.start = 10), 
          firstSelection = "top", depthCol = "nFrags", varFeatures = 25000, 
          dimsToUse = 1:30, LSIMethod = 2, scaleDims = TRUE, corCutOff = 0.75, 
          binarize = TRUE, outlierQuantiles = c(0.02, 0.98), filterBias = TRUE, 
          sampleCellsPre = 10000, projectCellsPre = FALSE, sampleCellsFinal = NULL, 
          selectionMethod = "var", scaleTo = 10000, totalFeatures = 5e+05, 
          filterQuantile = 0.995, excludeChr = c(), saveIterations = TRUE, 
          UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", 
                            verbose = FALSE, fast_sgd = TRUE), nPlot = 10000, outDir = getOutputDirectory(ArchRProj), 
          threads = getArchRThreads(), seed = 1, verbose = TRUE, force = FALSE, 
          logFile = createLogFile("addIterativeLSI")) 
{
  ArchR:::.startLogging(logFile = logFile)
  
  if (verbose) 
    message("Checking Inputs...")
  if (varFeatures < 1000) {
    stop("Please provide more than 1000 varFeatures!")
  }
  ArchR:::.requirePackage("Matrix", source = "cran")
  tstart <- Sys.time()
  if (!is.null(ArchRProj@reducedDims[[name]])) {
    if (!force) {
      stop("Error name in reducedDims Already Exists! Set force = TRUE or pick a different name!")
    }
  }
  set.seed(seed)
  outDir <- file.path(outDir, name)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  cellNames <- rownames(getCellColData(ArchRProj))
  if (!is.null(sampleCellsPre)) {
    if (length(cellNames) < sampleCellsPre) {
      sampleCellsPre <- NULL
    }
  }
  if (!is.null(sampleCellsFinal)) {
    if (length(cellNames) < sampleCellsFinal) {
      sampleCellsFinal <- NULL
    }
  }
  ArrowFiles <- getSampleColData(ArchRProj)[, "ArrowFiles"]
  if (tolower(useMatrix) == "tilematrix") {
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x) {
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if (length(unique(tileSizes)) != 1) {
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }
  else if (tolower(useMatrix) == "peakmatrix") {
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }
  else {
    tileSize <- NA
  }
  units <- unique(unlist(lapply(ArrowFiles, function(x) h5read(x, 
                                                               paste0(useMatrix, "/Info/Units")))))
  if (length(units) != 1) {
    stop("Units of matrices are not identical!")
  }
  if (grepl("log", units, ignore.case = TRUE)) {
    stop("Cannot use log transformed values for iterativeLSI!")
  }
  tstart <- Sys.time()
  chrToRun <- ArchR:::.availableSeqnames(ArrowFiles, subGroup = useMatrix)
  if (tolower(firstSelection) == "top") {
    if (!binarize) {
      stop("Please binarize data if using top selection for first iteration! Set binarize = TRUE!")
    }
    if (useMatrix == "TileMatrix") {
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, 
                              useMatrix = useMatrix, seqnames = chrToRun, 
                              addInfo = FALSE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }
    else {
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, 
                              useMatrix = useMatrix, seqnames = chrToRun, 
                              addInfo = TRUE)
    }
    if (length(excludeChr) > 0) {
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% 
                                                 excludeChr), , drop = FALSE]
    }
    nFeature <- varFeatures[1]
    rmTop <- floor((1 - filterQuantile) * totalFeatures)
    topIdx <- head(order(totalAcc$rowSums, decreasing = TRUE), 
                   nFeature + rmTop)[-seq_len(rmTop)]
    topFeatures <- totalAcc[sort(topIdx), ]
    gc()
  }
  else if (tolower(firstSelection) %in% c("var", "variable")) {
    if (binarize) {
      stop("Please do not binarize data if using variable selection for first iteration! Set binarize = FALSE!")
    }
    if (units %in% "BinarizedCounts") {
      stop("Cannot do variable selection with BinarizedCounts. Set firstSelection = Top!")
    }
    if (useMatrix == "TileMatrix") {
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, 
                              useMatrix = useMatrix, seqnames = chrToRun, 
                              useLog2 = TRUE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }
    else {
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, 
                              useMatrix = useMatrix, seqnames = chrToRun, 
                              useLog2 = TRUE)
    }
    if (length(excludeChr) > 0) {
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% 
                                                 excludeChr), , drop = FALSE]
    }
    nFeature <- varFeatures[1]
    if (nFeature > 0.5 * nrow(totalAcc)) {
      stop("nFeature for variable selection must be at leat 1/2 the total features!")
    }
    topIdx <- head(order(totalAcc$combinedVars, decreasing = TRUE), 
                   nFeature)
    topFeatures <- totalAcc[sort(topIdx), ]
    gc()
  }
  else {
    stop("firstSelect method must be Top or Var/Variable!")
  }
  cellDepth <- tryCatch({
    df <- getCellColData(ArchRProj = ArchRProj, select = depthCol)
    v <- df[, 1]
    names(v) <- rownames(df)
    v
  }, error = function(e) {
    tryCatch({
      ArchR:::.getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, 
                  seqnames = chrToRun)
    }, error = function(y) {
      stop("Could not determine depth from depthCol or colSums!")
    })
  })
  cellDepth <- log10(cellDepth + 1)
  j <- 1
  if (!is.null(clusterParams$sampleCells)) {
    if (!is.na(clusterParams$sampleCells[j])) {
      sampleJ <- clusterParams$sampleCells[j]
    }
    else if (!is.na(clusterParams$sampleCells[1])) {
      sampleJ <- clusterParams$sampleCells[1]
    }
    else {
      sampleJ <- sampleCellsPre
    }
  }
  else {
    sampleJ <- sampleCellsPre
  }
  outLSI <- ArchR:::.LSIPartialMatrix(ArrowFiles = ArrowFiles, featureDF = topFeatures, 
                              cellNames = cellNames, cellDepth = cellDepth, useMatrix = useMatrix, 
                              sampleNames = getCellColData(ArchRProj)$Sample, LSIMethod = LSIMethod, 
                              scaleTo = scaleTo, dimsToUse = dimsToUse, binarize = binarize, 
                              outlierQuantiles = outlierQuantiles, sampleCells = if (j != 
                                                                                     iterations) 
                                sampleCellsPre
                              else sampleCellsFinal, projectAll = j == iterations | 
                                projectCellsPre | sampleJ > sampleCellsPre, threads = threads, 
                              useIndex = FALSE, seed = seed, tstart = tstart, verbose = verbose, 
                              logFile = logFile)
  outLSI$scaleDims <- scaleDims
  outLSI$useMatrix <- useMatrix
  outLSI$tileSize <- tileSize
  gc()
  if (iterations == 1) {
    ArchRProj@reducedDims[[name]] <- outLSI
    return(ArchRProj)
  }
  clusterDF <- ArchR:::.LSICluster(outLSI = outLSI, filterBias = filterBias, 
                           cellNames = cellNames, cellDepth = cellDepth, dimsToUse = dimsToUse, 
                           scaleDims = scaleDims, corCutOff = corCutOff, clusterParams = clusterParams, 
                           j = j, verbose = verbose, tstart = tstart, logFile = logFile)
  clusters <- clusterDF$clusters
  nClust <- length(unique(clusters))
  if (saveIterations) {
    my_.saveIteration(outLSI = outLSI, clusters = clusters, 
                   scaleDims = scaleDims, dimsToUse = dimsToUse, corCutOff = corCutOff, 
                   outDir = outDir, nPlot = nPlot, UMAPParams = UMAPParams, 
                   ArchRProj = ArchRProj, j = j, threads = threads, 
                   logFile = logFile)
  }
  variableFeatures <- topFeatures
  while (j < iterations) {
    j <- j + 1
    variableFeatures <- ArchR:::.identifyVarFeatures(outLSI = outLSI, 
                                             clusterDF = clusterDF, ArrowFiles = ArrowFiles, 
                                             useMatrix = useMatrix, prevFeatures = variableFeatures, 
                                             scaleTo = scaleTo, totalAcc = totalAcc, totalFeatures = totalFeatures, 
                                             firstSelection = firstSelection, selectionMethod = selectionMethod, 
                                             varFeatures = varFeatures, tstart = tstart, threads = threads, 
                                             verbose = verbose, logFile = logFile)
    if (!is.null(clusterParams$sampleCells)) {
      if (!is.na(clusterParams$sampleCells[j])) {
        sampleJ <- clusterParams$sampleCells[j]
      }
      else if (!is.na(clusterParams$sampleCells[1])) {
        sampleJ <- clusterParams$sampleCells[1]
      }
      else {
        sampleJ <- sampleCellsPre
      }
    }
    else {
      sampleJ <- sampleCellsPre
    }
    outLSI <- ArchR:::.LSIPartialMatrix(ArrowFiles = ArrowFiles, 
                                featureDF = variableFeatures, useMatrix = useMatrix, 
                                cellNames = cellNames, cellDepth = cellDepth, sampleNames = getCellColData(ArchRProj)$Sample, 
                                LSIMethod = LSIMethod, scaleTo = scaleTo, dimsToUse = dimsToUse, 
                                binarize = binarize, outlierQuantiles = outlierQuantiles, 
                                sampleCells = if (j != iterations) 
                                  sampleCellsPre
                                else sampleCellsFinal, projectAll = j == iterations | 
                                  projectCellsPre | sampleJ > sampleCellsPre, 
                                threads = threads, useIndex = FALSE, seed = seed, 
                                tstart = tstart, verbose = verbose, logFile = logFile)
    outLSI$scaleDims <- scaleDims
    outLSI$useMatrix <- useMatrix
    outLSI$tileSize <- tileSize
    if (j != iterations) {
      clusterDF <- ArchR:::.LSICluster(outLSI = outLSI, dimsToUse = dimsToUse, 
                               scaleDims = scaleDims, corCutOff = corCutOff, 
                               filterBias = filterBias, cellNames = cellNames, 
                               cellDepth = cellDepth, j = j, clusterParams = clusterParams, 
                               verbose = verbose, tstart = tstart, logFile = logFile)
      clusters <- clusterDF$clusters
      nClust <- length(unique(clusters))
      if (saveIterations) {
        ArchR:::.saveIteration(outLSI = outLSI, clusters = clusters, 
                       scaleDims = scaleDims, dimsToUse = dimsToUse, 
                       corCutOff = corCutOff, outDir = outDir, nPlot = nPlot, 
                       UMAPParams = UMAPParams, ArchRProj = ArchRProj, 
                       j = j, threads = threads, logFile = logFile)
      }
    }
    print(dev.list())
    gc()
  }
  ArchRProj@reducedDims[[name]] <- outLSI
  return(ArchRProj)
}

my_.saveIteration <- function (outLSI = NULL, clusters = NULL, nPlot = NULL, UMAPParams = NULL, 
          ArchRProj = NULL, scaleDims = NULL, corCutOff = NULL, dimsToUse = NULL, 
          j = NULL, threads = NULL, outDir = NULL, logFile = NULL) 
{
  errorList <- append(args, mget(names(formals()), sys.frame(sys.nframe())))
  o <- tryCatch({
    ArchR:::.requirePackage("uwot", source = "cran")
    if (scaleDims) {
      dimsPF <- dimsToUse[which(outLSI$corToDepth$scaled[dimsToUse] <= 
                                  corCutOff)]
    }
    else {
      dimsPF <- dimsToUse[which(outLSI$corToDepth$none[dimsToUse] <= 
                                  corCutOff)]
    }
    if (nrow(outLSI[[1]]) > nPlot) {
      saveIdx <- sample(seq_len(nrow(outLSI[[1]])), nPlot)
    }
    else {
      saveIdx <- seq_len(nrow(outLSI[[1]]))
    }
    UMAPParams <- ArchR:::.mergeParams(UMAPParams, list(n_neighbors = 40, 
                                                min_dist = 0.4, metric = "cosine", verbose = FALSE, 
                                                fast_sgd = TRUE))
    if (scaleDims) {
      UMAPParams$X <- ArchR:::.scaleDims((outLSI[[1]][saveIdx, 
                                              , drop = FALSE])[, dimsPF, drop = FALSE])
    }
    else {
      UMAPParams$X <- (outLSI[[1]][saveIdx, , drop = FALSE])[, 
                                                             dimsPF, drop = FALSE]
    }
    UMAPParams$ret_nn <- FALSE
    UMAPParams$ret_model <- FALSE
    UMAPParams$n_threads <- floor(threads/2)
    uwotUmap <- do.call(uwot::umap, UMAPParams)
    cat("flag 1: ", dev.cur(), "\n")
    pdf(file.path(outDir, paste0("Save-LSI-Iteration-", 
                                 j, ".pdf")), width = 6, height = 6)
    cat("flag 2: ", dev.cur(), "\n")
    p1 <- my_ggPoint(uwotUmap[, 1], uwotUmap[, 2], getCellColData(ArchRProj, 
                                                               select = "Sample")[rownames(outLSI[[1]])[saveIdx], 
                                                               ], size = 0.5, title = paste0("SampleName (nCells Plot = ", 
                                                                                             nrow(UMAPParams$X), ")"), rastr = TRUE)
    p2 <- my_ggPoint(uwotUmap[, 1], uwotUmap[, 2], clusters[saveIdx], 
                  size = 0.5, title = paste0("Clusters (nCells Plot = ", 
                                             nrow(UMAPParams$X), ")"), rastr = TRUE)
    cat("flag 3: ", dev.cur(), "\n")
    p1 <- p1 + xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p2 <- p2 + xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
    cat("flag 4: ", dev.cur(), "\n")
    ArchR:::.fixPlotSize(p1, plotWidth = 6, plotHeight = 6)
    grid::grid.newpage()
    ArchR:::.fixPlotSize(p2, plotWidth = 6, plotHeight = 6)
    dev.off()
    outj <- SimpleList(LSI = outLSI, clusters = clusters, 
                       uwotUmap = uwotUmap)
    ArchR:::.safeSaveRDS(outj, file.path(outDir, paste0("Save-LSI-Iteration-", 
                                                j, ".rds")))
    rm(UMAPParams, uwotUmap)
    gc()
  }, error = function(e) {
    print(e)
  })
  return(0)
}

my_ggPoint <- function (x = NULL, y = NULL, color = NULL, discrete = TRUE, 
          discreteSet = "stallion", continuousSet = "solarExtra", 
          labelMeans = TRUE, pal = NULL, defaultColor = "lightGrey", 
          highlightPoints = NULL, colorDensity = FALSE, size = 1, 
          xlim = NULL, ylim = NULL, extend = 0.05, xlabel = "x", ylabel = "y", 
          title = "", randomize = FALSE, seed = 1, colorTitle = NULL, 
          colorOrder = NULL, colorLimits = NULL, alpha = 1, baseSize = 10, 
          legendSize = 3, ratioYX = 1, labelAsFactors = TRUE, fgColor = "black", 
          bgColor = "white", bgWidth = 1, labelSize = 3, addFit = NULL, 
          rastr = FALSE, dpi = 300, ...) 
{
  stopifnot(length(y) == length(x))
  if (length(x) < 5) {
    stop("x must be at least length 5 to plot!")
  }
  if (randomize) {
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  }
  else {
    idx <- seq_along(x)
  }
  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  if (length(include) != length(x)) {
    message("Some values are not finite! Excluding these points!")
    df <- df[include, ]
    x <- x[include]
    y <- y[include]
    if (!is.null(color)) {
      color <- color[include]
    }
  }
  if (is.null(xlim)) {
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  if (is.null(ylim)) {
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  cat("in 1: ", dev.cur(), "\n")
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  ArchR:::.requirePackage("ggplot2", source = "cran")
  cat("in 2: ", dev.cur(), "\n")
  if (is.null(color) & !colorDensity) {
    p <- ggplot(df[idx, ], aes(x = x, y = y)) + coord_equal(ratio = ratioXY, 
                                                            xlim = xlim, ylim = ylim, expand = F) + xlab(xlabel) + 
      ylab(ylabel) + ggtitle(title) + theme_ArchR(baseSize = baseSize)
    if (rastr) {
      p <- p + ArchR:::.geom_point_rast2(size = size, raster.dpi = dpi, 
                                 alpha = alpha, color = defaultColor)
    }
    else {
      p <- p + geom_point(size = size, alpha = alpha, 
                          color = defaultColor)
    }
  }
  else {
    if (colorDensity) {
      discrete <- FALSE
      df <- ArchR:::.getDensity(x, y, n = 100, sample = NULL)
      df <- df[order(df$density), , drop = FALSE]
      df$color <- df$density
      if (is.null(colorTitle)) {
        colorTitle <- "density"
      }
    }
    else if (discrete) {
      if (!is.null(highlightPoints)) {
        if (length(highlightPoints) < length(color)) {
          color[-highlightPoints] <- "Non.Highlighted"
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      color <- paste0(color)
      if (!is.null(colorOrder)) {
        if (!all(color %in% colorOrder)) {
          stop("Not all colors are in colorOrder!")
        }
      }
      else {
        colorOrder <- gtools::mixedsort(unique(color))
      }
      if (is.null(colorTitle)) {
        colorTitle <- "color"
      }
      stopifnot(length(color) == nrow(df))
      df$color <- factor(color, levels = colorOrder)
      if (labelAsFactors) {
        df$color <- factor(x = paste0(paste0(match(paste0(df$color), 
                                                   paste0(levels(df$color)))), "-", paste0(df$color)), 
                           levels = paste0(seq_along(levels(df$color)), 
                                           "-", levels(df$color)))
        if (!is.null(pal)) {
          names(pal) <- paste0(levels(df$color))[match(names(pal), 
                                                       colorOrder)]
        }
        colorOrder <- paste0(levels(df$color))
      }
    }
    else {
      stopifnot(length(color) == nrow(df))
      if (!is.null(highlightPoints)) {
        if (length(highlightPoints) < length(color)) {
          color[-highlightPoints] <- NA
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      if (!is.null(colorLimits)) {
        color[color < min(colorLimits)] <- min(colorLimits)
        color[color > max(colorLimits)] <- max(colorLimits)
      }
      df$color <- color
    }
    p <- ggplot(df[idx, ], aes(x = x, y = y, color = color)) + 
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, 
                  expand = FALSE) + xlab(xlabel) + ylab(ylabel) + 
      ggtitle(title) + theme_ArchR(baseSize = baseSize) + 
      theme(legend.direction = "horizontal", legend.box.background = element_rect(color = NA)) + 
      labs(color = colorTitle)
    if (rastr) {
      p <- p + ArchR:::.geom_point_rast2(size = size, raster.dpi = dpi, 
                                 alpha = alpha, raster.width = min(par("fin")), 
                                 raster.height = (ratioYX * min(par("fin"))))
    }
    else {
      p <- p + geom_point(size = size, alpha = alpha)
    }
    if (discrete) {
      if (!is.null(pal)) {
        p <- p + scale_color_manual(values = pal)
      }
      else {
        pal <- paletteDiscrete(set = discreteSet, values = colorOrder)
        if (!is.null(highlightPoints)) {
          pal[grep("Non.Highlighted", names(pal))] <- "lightgrey"
        }
        p <- p + scale_color_manual(values = pal) + 
          guides(color = guide_legend(override.aes = list(size = legendSize, 
                                                          shape = 15)))
      }
      if (labelMeans) {
        dfMean <- split(df, df$color) %>% lapply(., 
                                                 function(x) {
                                                   data.frame(x = median(x[, 1]), y = median(x[, 
                                                                                               2]), color = x[1, 3])
                                                 }) %>% Reduce("rbind", .)
        if (labelAsFactors) {
          dfMean$label <- stringr::str_split(paste0(seq_len(nrow(dfMean))), 
                                             pattern = "\\-", simplify = TRUE)[, 1]
        }
        else {
          dfMean$label <- dfMean$color
        }
        dfMean$text <- stringr::str_split(dfMean$color, 
                                          pattern = "-", simplify = TRUE)[, 1]
        theta <- seq(pi/8, 2 * pi, length.out = 16)
        xo <- bgWidth * diff(range(df$x))/300
        yo <- bgWidth * diff(range(df$y))/300
        for (i in theta) {
          p <- p + geom_text(data = dfMean, aes_q(x = bquote(x + 
                                                               .(cos(i) * xo)), y = bquote(y + .(sin(i) * 
                                                                                                   yo)), label = ~text), size = labelSize, 
                             color = bgColor)
        }
        if (is.null(fgColor)) {
          p <- p + geom_text(data = dfMean, aes(x = x, 
                                                y = y, color = color, label = label), size = labelSize, 
                             show.legend = FALSE)
        }
        else {
          p <- p + geom_text(data = dfMean, aes(x = x, 
                                                y = y, label = label), color = fgColor, 
                             size = labelSize, show.legend = FALSE)
        }
      }
    }
    else {
      if (!is.null(pal)) {
        if (!is.null(colorLimits)) {
          p <- p + scale_colour_gradientn(colors = pal, 
                                          limits = colorLimits, na.value = "lightgrey")
        }
        else {
          p <- p + scale_colour_gradientn(colors = pal, 
                                          na.value = "lightgrey")
        }
      }
      else {
        if (!is.null(colorLimits)) {
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), 
                                          limits = colorLimits, na.value = "lightgrey")
        }
        else {
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), 
                                          na.value = "lightgrey")
        }
      }
    }
  }
  cat("in 3: ", dev.cur(), "\n")
  if (!is.null(addFit)) {
    p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, 
                         color = "black") + ggtitle(paste0(title, "\nPearson = ", 
                                                           round(cor(df$x, df$y), 3), "\nSpearman = ", round(cor(df$x, 
                                                                                                                 df$y, method = "spearman"), 3)))
  }
  cat("in 4: ", dev.cur(), "\n")
  p <- p + theme(legend.position = "bottom", legend.key = element_rect(size = 2))
  cat("in 5: ", dev.cur(), "\n")
  if (!is.null(ratioYX)) {
    attr(p, "ratioYX") <- ratioYX
  }
  return(p)
}







# 3 use peak called by ArchR iterative pipeline
# 3.1 peak calling
proj <- addGroupCoverages(proj, groupBy = "predictedIdent", minReplicates = 2)

proj <- addReproduciblePeakSet(proj, groupBy = "predictedIdent", force = T, peakMethod = "Macs2", 
                               pathToMacs2 = "/home/chenzy/miniconda/envs/macs2/bin/macs2")



# Calculate correlation between cicero and archR gene score
dim(unnorm_ga)
gs_mat <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
dim(gs_mat)
gene_name <- elementMetadata(gs_mat)$name
gs_mat <- assay(gs_mat)
rownames(gs_mat) <- gene_name

genes_use <- intersect(rownames(gs_mat), rownames(unnorm_ga))
cell_order_list <- colnames(unnorm_ga)

unnorm_ga <- unnorm_ga[genes_use, ]
gs_mat <- gs_mat[genes_use, cell_order_list]

corr <- sapply(1:length(colnames(unnorm_ga)), FUN = function(i) {cor(unnorm_ga[, i], gs_mat[, i])})











# try to implement self peak adding
addMyPeakset <- function(proj, groupBy = "Clusters",
                         peakMethod = "Macs2",
                         reproducibility = "2",
                         peaksPerCell = 500,
                         maxPeaks = 150000,
                         minCells = 25,
                         excludeChr = c("chrM","chrY"),
                         pathToMacs2 = if(tolower(peakMethod)=="macs2") findMacs2() else NULL,
                         genomeSize = NULL, 
                         shift = -75, 
                         extsize = 150, 
                         method = if(tolower(peakMethod)=="macs2") "q" else "p", #P-Method for Tiles Results Better Agree w/ Macs2
                         cutOff = 0.1, 
                         additionalParams = "--nomodel --nolambda",
                         extendSummits = 250,
                         promoterRegion = c(2000, 100),
                         genomeAnnotation = getGenomeAnnotation(ArchRProj),
                         geneAnnotation = getGeneAnnotation(ArchRProj),
                         plot = TRUE,
                         threads = getArchRThreads(),
                         parallelParam = NULL,
                         force = FALSE,
                         verbose = TRUE,
                         logFile = createLogFile("addReproduciblePeakSet")){
  coverageMetadata <- .getCoverageMetadata(ArchRProj = ArchRProj, groupBy = groupBy, minCells = minCells)
  coverageParams <- .getCoverageParams(ArchRProj = ArchRProj, groupBy = groupBy)
  
  tableGroups <- table(getCellColData(ArchRProj, groupBy, drop = TRUE))
  groupSummary <- lapply(seq_along(coverageParams$cellGroups), function(y){
    x <- coverageParams$cellGroups[[y]]
    uniq <- unique(unlist(x))
    n <- lapply(x, length) %>% unlist %>% sum
    nmin <- lapply(x, length) %>% unlist %>% min
    nmax <- lapply(x, length) %>% unlist %>% max
    data.frame(
      Group=names(coverageParams$cellGroups)[y], 
      nCells=tableGroups[names(coverageParams$cellGroups)[y]], 
      nCellsUsed=length(uniq), 
      nReplicates=length(x), 
      nMin=nmin, 
      nMax=nmax, 
      maxPeaks = min(maxPeaks, length(uniq) * peaksPerCell)
    )
  }) %>% Reduce("rbind",.)
  
  #####################################################
  # Create Output Directory
  #####################################################
  outDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls")
  outSubDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls", "ReplicateCalls")
  outBedDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls", "InsertionBeds")
  dir.create(outDir, showWarnings = FALSE)
  dir.create(outSubDir, showWarnings = FALSE)
  dir.create(outBedDir, showWarnings = FALSE)
  
  #####################################################
  # Genome Size Presets
  #####################################################
  if(is.null(genomeSize)){
    if(grepl("hg19|hg38", getGenome(ArchRProj), ignore.case = TRUE)){
      genomeSize <- 2.7e9
    }else if(grepl("mm9|mm10", getGenome(ArchRProj), ignore.case = TRUE)){
      genomeSize <- 1.87e9
    }
  }
  
  #####################################################
  # Arguments for Peak Calling
  #####################################################
  coverageFiles <- coverageMetadata$File
  names(coverageFiles) <- coverageMetadata$Name
  args <- list()
  args$X <- seq_len(nrow(coverageMetadata))
  args$FUN <- .callSummitsOnCoverages
  args$coverageFiles <- coverageFiles
  args$outFiles <- file.path(outSubDir, paste0(make.names(coverageMetadata$Name),"-summits.rds"))
  args$bedDir <- outBedDir
  args$excludeChr <- excludeChr
  args$peakParams <- list(
    pathToMacs2 = pathToMacs2,
    genomeSize = genomeSize, 
    shift = shift, 
    extsize = extsize, 
    cutOff = cutOff, 
    method = method,
    additionalParams = additionalParams
  )
  args$parallelParam <- parallelParam
  args$threads <- threads
  # args$logFile <- logFile
  args$registryDir <- file.path(outDir, "batchRegistry")
  
  #back lapply
  outSummmits <- unlist(.batchlapply(args))
  
  #Summarize Output
  outSummitList <- split(outSummmits, coverageMetadata$Group)
  summitNamesList <- split(coverageMetadata$Name, coverageMetadata$Group)
  
  outSummmits <- unlist(.batchlapply(args))
  
}

