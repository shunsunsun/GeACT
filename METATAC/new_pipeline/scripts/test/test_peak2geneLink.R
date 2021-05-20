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

proj_epi <- loadArchRProject("ArchR_epi", showLogo = F)


ArchRProj = proj_epi
reducedDims = "peakLSI"
useMatrix = "GeneIntegrationMatrix"
dimsToUse = 1:30
scaleDims = NULL
corCutOff = 0.75
cellsToUse = NULL
k = 100
knnIteration = 500
overlapCutoff = 0.8
maxDist = 25000
scaleTo = 10^4
log2Norm = TRUE
predictionCutoff = 0.4
addEmpiricalPval = FALSE
seed = 1
threads = max(floor(getArchRThreads() / 2), 1)
verbose = TRUE
logFile = NULL

addPeak2GeneLinks <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  addEmpiricalPval = FALSE,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = NULL
){
  tstart <- Sys.time()
  AvailableMatrices <- getAvailableMatrices(ArchRProj)
  
  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }
  
  if(useMatrix %ni% AvailableMatrices){
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  
  
  dfAll <- ArchR:::.safelapply(seq_along(ArrowFiles), function(x){
    cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames")))
    pSx <- tryCatch({
      h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    }, error = function(e){
      if(getArchRVerbose()) message("No predictionScore found. Continuing without predictionScore!")
      rep(9999999, length(cNx))
    })
    DataFrame(
      cellNames = cNx,
      predictionScore = pSx
    )
  }, threads = threads) %>% Reduce("rbind", .)
  
  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]
  
  set.seed(seed)
  
  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  
  #Gene Info
  geneSet <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  
  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }
  
  #Subsample (select subset of cells to compute k nearest neighbors)
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
  
  #KNN Matrix
  knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx,], k = k)
  
  #Determin Overlap (omit cells that have too many overlaps in their nearest neighbors, that is to say, keep the dissimilarity of selected cells)
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  
  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  #Check Chromosomes for peaks and genes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)
  
  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)
  
  #Group Matrix RNA
  groupMatRNA <- ArchR:::.getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  rawMatRNA <- groupMatRNA
  
  #Group Matrix ATAC
  groupMatATAC <- ArchR:::.getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  rawMatATAC <- groupMatATAC
  
  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo
  
  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }
  
  names(geneStart) <- NULL
  
  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj
  
  names(peakSet) <- NULL
  
  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
    rowRanges = peakSet
  )
  metadata(seATAC)$KNNList <- knnObj
  
  rm(groupMatRNA, groupMatATAC)
  gc()
  
  #Overlaps between gene and peak granges (by indexes)
  o <- DataFrame(
    findOverlaps(
      ArchR:::.suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), 
      ignore.strand = TRUE
    )
  )
  
  #Get Distance from Fixed point A B, ???? what is distance used for?
  o$distance <- distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] )
  colnames(o) <- c("B", "A", "distance")
  
  #Null Correlations
  if(addEmpiricalPval){
    # "Computing Background Correlations"
    nullCor <- ArchR:::.getNullCorrelations(seATAC, seRNA, o, 1000)
  }
  
  # "Computing Correlations"
  o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 1e-17, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")  
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart
  
  if(addEmpiricalPval){
    out$EmpPval <- 2*pnorm(-abs(((out$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }
  
  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  ArchR:::.safeSaveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  ArchR:::.safeSaveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA
  
  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out
  
  ArchRProj
  
}