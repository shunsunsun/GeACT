# peak overlapping with conserved elements (ATAC)
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("ggplot2")
library("cowplot")

ovlpStat <- read.table(file = "04-open_chromatin/merged/ovlpStat.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
colnames(ovlpStat) <- c("q", "m", "n", "k")
ovlpStat$a <- ovlpStat$m + ovlpStat$n
ovlpStat$bgRatio <- ovlpStat$m / ovlpStat$a
ovlpStat$pkRatio <- ovlpStat$q / ovlpStat$k
ovlpStat$region <- rownames(ovlpStat)

cat("Ratio of the peak regions in non-CDS region:", ovlpStat["non-CDS", "k"] / ovlpStat["All", "k"], "\n")

apply(ovlpStat[, 1:7], 1, function(x) {
  phyper(q = x[1], m = x[2], n = x[3], k = x[4], lower.tail = F)
})

stat_DF <- melt(ovlpStat[, c("region", "bgRatio", "pkRatio")], id.vars = "region", variable.name = "type", value.name = "consRatio")
levels(stat_DF$type)
levels(stat_DF$type) <- c("Genome", "All peaks")
levels(stat_DF$type)
stat_DF$consPercentage <- stat_DF$consRatio * 100
#
stat_DF <- subset(stat_DF, region == "All")
#
pdf(file = "peakConsRatio.pdf", width = 5, height = 5)
ymax <- max(stat_DF$consPercentage) * 1.2
sigBar <- data.frame(type = rep(levels(stat_DF$type), each = 2), consPercentage = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(stat_DF, aes(x = type, y = consPercentage, fill = type)) + geom_col(show.legend = F) + 
  #facet_grid(. ~ region) + 
  geom_path(data = sigBar, group = 1) + geom_text(x = 1.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(limits = c(0, 0.13 * 100), expand = c(0, 0)) + 
  scale_fill_manual(values = c("gray75", "cornflowerblue")) + 
  xlab(NULL) + ylab("Percentage of conserved regions") + 
  theme(aspect.ratio = 1)
dev.off()
