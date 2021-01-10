
do_cellFiltering <- function(do_plot = T) {
  filtering_DF <- data.frame(ABratio=cellStat$ABratio >= cutoff_DF$ABratio, cleanReads=cellStat$cleanReads >= cutoff_DF$cleanReads, 
                             mpRatio=cellStat$mpRatio >= cutoff_DF$mpRatio, 
                             nGene_l=cellStat$nGene > cutoff_DF$nGene_l, nGene_h=cellStat$nGene < cutoff_DF$nGene_h, 
                             nUMI_l=cellStat$nUMI > cutoff_DF$nUMI_l, nUMI_h=cellStat$nUMI < cutoff_DF$nUMI_h, 
                             mitoRatio=cellStat$mitoRatio < cutoff_DF$mitoRatio, ERCCratio=cellStat$ERCCratio < cutoff_DF$ERCCratio, 
                             avgCor=cellStat$avgCor >= cutoff_DF$avgCor, 
                             doublet= ! cellStat$doublet)
  rownames(filtering_DF) <- rownames(cellStat)
  filtering_MT <- apply(filtering_DF, 1, function(x) { y <- ifelse(x, 1, 0); return(y) })
  hc <- hclust(dist(t(filtering_MT)))
  filtering_MT_melted <- melt(data = filtering_MT)
  colnames(filtering_MT_melted) <- c("index", "cell", "value")
  filtering_MT_melted$cell <- factor(filtering_MT_melted$cell, levels = colnames(filtering_MT)[hc$order])
  # plot
  if(do_plot) {
    p <- ggplot(filtering_MT_melted, aes(x = cell, y = index, fill = as.factor(value))) + geom_tile(show.legend = F) + theme_bw() + 
      theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(size = 1, color = "black")) + 
      scale_y_discrete(limits = rev(levels(filtering_MT_melted$index)), position = "right", expand = c(0, 0))
    print(p)
  }
  # stat
  filtering_res <- apply(filtering_DF, 1, all)
  filtering_res_ratio <- sum(filtering_res) / length(filtering_res)
  cat(sum(filtering_res), "/" , length(filtering_res), "cells remained:", round(filtering_res_ratio * 100, 2), "%\n")
  
  out <- list(res = filtering_res, tb = filtering_MT_melted)
  return(out)
}

do_plotIndex <- function() {
  par(mfrow=c(4,2), mar = c(4.5, 4.2, 0.5, 0.5))
  hist(cellStat$ABratio, breaks = 40, xlab = "A/B ratio", main = NA)
  hist(cellStat$cleanReads/1e6, breaks = 40, xlab = "Clean reads (M)", main = NA)
  hist(cellStat$mpRatio, breaks = 40, xlab = "Mapping ratio", main = NA, xlim = c(0,1))
  hist(cellStat$nGene, breaks = 40, xlab = "Detected genes", main = NA)
  hist(cellStat$nUMI/1e3, breaks = 80, xlab = "Detected transcripts (k)", main = NA)
  hist(cellStat$mitoRatio, breaks = 40, xlab = "Mito. ratio", main = NA)
  hist(cellStat$ERCCratio, breaks = 40, xlab = "ERCC ratio", main = NA)
  hist(cellStat$avgCor, breaks = 40, xlab = "Pairwise correlation", main = NA)
  par(mfrow=c(1,1), mar = c(5,4,4,2) + 0.1)
}

do_show_ftd_stat <- function(cellStat) {
  out1 <- data.frame(total = length(cellStat$filter), fail = sum(! cellStat$filter), pass = sum(cellStat$filter), 
             ratio = sum(cellStat$filter) / length(cellStat$filter), row.names = "cell")
  print(out1)
  cat("\n")
  out2 <- data.frame(average = c(round(mean(cellStat[cellStat$filter, "cleanReads"]/1e6), 2), 
                         round(mean(cellStat[cellStat$filter, "nGene"]), 2), 
                         round(mean(cellStat[cellStat$filter, "nUMI"]), 2)), 
             median = c(round(median(cellStat[cellStat$filter, "cleanReads"]/1e6), 2), 
                        round(median(cellStat[cellStat$filter, "nGene"]), 2), 
                        round(median(cellStat[cellStat$filter, "nUMI"]), 2)), 
             row.names = c("cleanReads", "nGene", "nUMI"))
  print(out2)
}
