library(ggplot)
library(tidyverse)
library(cowplot)
library(GenomicRanges)
library("VennDiagram")

setwd("/data/Lab/otherwork/GeACT/ATAC/19-21w/02_small_intestine/")

peak_all <- "../../macs2/all_organ_peaks.narrowPeak"
peak_file1 <- "peak_calling/macs2_legacy/02_small_intestine_peaks.narrowPeak"
peak_file2 <- "peak_calling/macs2/organ_peaks.narrowPeak"


my_clip <- function(x, lb, ub){
  pmax(lb, pmin(x, ub))
}

# length distribution
plot_peak_dis <- function(peak_file, binwidth = 50, cutoff = 500){
  peak <- read_tsv(peak_file, col_names = F)
  
  print(table((peak$X3 - peak$X2) <= cutoff))
  
  # number of peaks
  print(sprintf("peaks count: %d", nrow(peak)))
  
  p1 <- ggplot(data = peak, mapping = aes(x = my_clip(X3 - X2, 0, 5000))) + geom_histogram(binwidth = binwidth)
  p2 <- ggplot(data = peak, mapping = aes(x = "peak length", y = my_clip(X3 - X2, 0, 5000))) + geom_boxplot()
  
  print(p1)
  
  return(list(p1, p2))
}

# overlap between two peakset
plot_intersection <- function(peak_file1, peak_file2){
  peak1 <- rtracklayer::import(peak_file1)
  peak2 <- rtracklayer::import(peak_file2)
  ipeak <- GenomicRanges::intersect(peak1, peak2) %>% GenomicRanges::reduce()
  
  peak1_list <- split(peak1, seqnames(peak1))
  peak2_list <- split(peak2, seqnames(peak2))
  ipeak_list <- split(ipeak, seqnames(ipeak))
  
  peak1_coverage <- sapply(peak1_list, function(gr){sum(width(gr))})
  peak2_coverage <- sapply(peak2_list, function(gr){sum(width(gr))})
  ipeak_coverage <- sapply(ipeak_list, function(gr){sum(width(gr))})
  
  print(peak1_coverage)
  print(peak2_coverage)
  print(ipeak_coverage)
}


figs0 <- plot_peak_dis(peak_all)
figs1 <- plot_peak_dis(peak_file1)
figs2 <- plot_peak_dis(peak_file2)

plot_intersection(peak_all, peak_file1)
plot_intersection(peak_file1, peak_file2)






