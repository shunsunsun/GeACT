#!/usr/bin/env r
suppressMessages({
  library(docopt)
})


# Argument parsing -------------------------------------------------------------
doc <- "Usage: merge_peakSet.R [-h] [-b] [-s] [-o OUTPUT] [-a ANNODIR] [PEAKS ...]

-h --help           show help message
-b --frombed        load peaks data from bed, default from rds
-s --simplemerge    simple merge, instead of iteratively merge
-o --output OUTPUT  output prefix **path** [default: ./union_stage]
-a --annodir ANNODIR   gene and genome annotation dir [default: /data/Lab/otherwork/GeACT/ATAC/database/annotation]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

suppressMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

frombed <- opt$frombed
simplemerge <- opt$simplemerge
annoDir <- opt$annodir
outprefix <- opt$output

# opt$PEAKS <- c("/data/Lab/otherwork/GeACT/ATAC/data/19-22w/01_stomach/results/peak_calling/normalPeaks.rds", 
#                "/data/Lab/otherwork/GeACT/ATAC/data/11-14w/01_stomach/results/peak_calling/normalPeaks.rds")
multi_peaks <- lapply(opt$PEAKS, function(peak_file){
  if(frombed){
    peaksDF <- read.table(peak_file, sep = "\t", quote = "", head = T)
    peaks <- GRanges(peaksDF$seqnames, IRanges(peaksDF$start + 1, peaksDF$end))
    mcols(peaks) <- peaksDF[, !(colnames(peaksDF) %in% c("seqnames", "start", "end", "strand", "width"))]
  }
  else
    peaks <- readRDS(peak_file)
  peaks$quantileScore <- round(trunc(rank(peaks$score))/length(peaks$score), 3)
  return(peaks)
}) %>% GRangesList

merge_peaks <- unlist(multi_peaks)

if (simplemerge){
  merge_peaks <- sort(sortSeqlevels(merge_peaks))
  merge_peaks <- reduce(merge_peaks)
  
  # load genome and gene annotation
  geneAnno <- readRDS(paste(annoDir, "geneAnnotation.rds", sep = "/"))
  genomeAnno <- readRDS(paste(annoDir, "genomeAnnotation.rds", sep = "/"))
  
  BSgenome <- eval(parse(text = genomeAnno$genome))
  BSgenome <- validBSgenome(BSgenome)
  promoterRegion <- c(2000, 100)
  
  merge_peaks <- ArchR:::.fastAnnoPeaks(merge_peaks, BSgenome = BSgenome, geneAnnotation = geneAnno, promoterRegion = promoterRegion)
  mcols(merge_peaks)$N <- NULL
} else{
  merge_peaks <- nonOverlappingGR(merge_peaks, by = "quantileScore", decreasing = T)
}
  

if (simplemerge){
  saveRDS(merge_peaks, file = paste0(outprefix, "_simplemerge.rds"))
  
  peaksDF <- data.frame(chr=seqnames(merge_peaks), start=start(merge_peaks) - 1, end=end(merge_peaks))
  peaksDF <- cbind(peaksDF, mcols(merge_peaks))
  write.table(peaksDF, file = paste0(outprefix, "_simplemerge.bed"), col.names = T, quote = F, row.names = F, sep = "\t")
} else {
  saveRDS(merge_peaks, file = paste0(outprefix, ".rds"))
  
  peaksDF <- data.frame(chr=seqnames(merge_peaks), start=start(merge_peaks) - 1, end=end(merge_peaks))
  peaksDF <- cbind(peaksDF, mcols(merge_peaks))
  write.table(peaksDF, file = paste0(outprefix, ".bed"), col.names = T, quote = F, row.names = F, sep = "\t")
}
