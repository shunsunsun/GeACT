#!/usr/bin/env r
suppressMessages({
  library(docopt)
})


# Argument parsing -------------------------------------------------------------
doc <- "Usage: construct_peakSet.R [-h] [--inputdir INDIR] [--output OUTPUT] [--source SOURCE] [--annodir ANNODIR] [--fromSummits] [--extendSummits EXTEND] [--maxPeaks MAXPEAKS]

-h --help           show help message
--inputdir INDIR    input directory containing summits or narrowPeaks file [default: /data/Lab/otherwork/GeACT/ATAC/data/19-22w/01_stomach/results/peak_calling/macs2]
--output OUTPUT     output prefix normalPeaks or summitPeaks
--source SOURCE     source name used for later merging [default: default]
--annodir ANNODIR   gene and genome annotation dir [default: /data/Lab/otherwork/GeACT/ATAC/database/annotation]
--fromSummits       construct peakSet by summits, default by narrowPeak
--extendSummits EXTEND length of extending summits, used with fromSummits [default: 250]
--maxPeaks MAXPEAKS max number of peaks to keep [default: 150000]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

inDir <- opt$inputdir
if(is.null(opt$output)) outputPrefix <- paste(inDir, "normalPeaks", sep = "/") else 
  outputPrefix <- opt$output

extendSummits <- as.integer(opt$extendSummits)
fromSummits <- opt$fromSummits
maxPeaks <- as.integer(opt$maxPeaks)
annoDir <- opt$annodir
sourceName <- opt$source

suppressMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# load genome and gene annotation
geneAnno <- readRDS(paste(annoDir, "geneAnnotation.rds", sep = "/"))
genomeAnno <- readRDS(paste(annoDir, "genomeAnnotation.rds", sep = "/"))

if(fromSummits){
  summitFile <- paste(inDir, "organ_summits.bed", sep = "/")
  summits <- rtracklayer::import(summitFile)
  extSummits <- resize(summits, extendSummits * 2 + 1, "center")
  
  extSummits <- subsetByOverlaps(extSummits, genomeAnno$blacklist, invert = TRUE)
  peaks <- nonOverlappingGR(extSummits, by = "score", decreasing = TRUE)
} else{
  narrowPeaksFile <- paste(inDir, "organ_peaks.narrowPeak", sep = "/")
  peaks <- rtracklayer::import(narrowPeaksFile)
  peaks <- nonOverlappingGR(peaks, by = "score", decreasing = TRUE)
}

BSgenome <- eval(parse(text = genomeAnno$genome))
BSgenome <- validBSgenome(BSgenome)
promoterRegion <- c(2000, 100)

peaks <- sort(sortSeqlevels(peaks))
peaks <- subsetByOverlaps(peaks, genomeAnno$chromSizes, type = "within")
peaks <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnno, promoterRegion = promoterRegion)

peaks <- peaks[which(mcols(peaks)$N < 0.001)] #Remove N Containing Peaks
peaks <- peaks[order(peaks$score, decreasing = TRUE)]
peaks <- head(peaks, maxPeaks)
mcols(peaks)$N <- NULL

peaks <- sort(peaks)
mcols(peaks)$name <- NULL
mcols(peaks)$sourceName <- sourceName

saveRDS(peaks, file = paste0(outputPrefix, ".rds"))

peaksDF <- data.frame(chr=seqnames(peaks), start=start(peaks) - 1, end=end(peaks))
peaksDF <- cbind(peaksDF, mcols(peaks))
write.table(peaksDF, file = paste0(outputPrefix, ".bed"), col.names = T, quote = F, row.names = F, sep = "\t")
