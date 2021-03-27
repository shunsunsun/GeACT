
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")
library("ggplot2")
library("cowplot")

OUT <- "04-open_chromatin/merged"

TF_enrich <- read.table(file = "04-open_chromatin/merged/module_TF_metatable_dedup.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(TF_enrich)

# MD91
mdid_case <- "MD91"
TF_enrich_sub <- subset(TF_enrich, MID == mdid_case & celltype == "Small intestine.Fibro-COL6A5")
TF_enrich_sub$q <- p.adjust(TF_enrich_sub$p)
TF_enrich_sub$highlight <- TF_enrich_sub$TF_module_avgCor > 0.15 & TF_enrich_sub$q < 1e-5
table(TF_enrich_sub$highlight)

pdf(file = paste0(OUT, "/", mdid_case, "_TF_avgCor.pdf"), width = 5, height = 4)

ggplot(TF_enrich_sub, aes(x = TF_module_avgCor, y = -log10(q))) + geom_point(aes(color = highlight), alpha = 0.6, show.legend = F) + 
  geom_text(data = subset(TF_enrich_sub, highlight), aes(label = TF), nudge_x = 0.01, nudge_y = 0.6) + 
  scale_color_manual(values = c("grey50", "tomato")) + 
  xlab("Average TF-gene expression correlation") + ylab(expression(paste(-Log[10], " (FDR)"))) + 
  coord_cartesian(ylim = c(0, 10), clip = "off") + 
  theme(aspect.ratio = 1)

dev.off()

# MD203
mdid_case <- "MD203"
TF_enrich_sub <- subset(TF_enrich, MID == mdid_case & celltype == "Small intestine.Fibro-COL6A5")
TF_enrich_sub$q <- p.adjust(TF_enrich_sub$p)
TF_enrich_sub$highlight <- TF_enrich_sub$TF_module_avgCor > 0.1 & TF_enrich_sub$q < 0.0005
table(TF_enrich_sub$highlight)

pdf(file = paste0(OUT, "/", mdid_case, "_TF_avgCor.pdf"), width = 5, height = 4)

ggplot(TF_enrich_sub, aes(x = TF_module_avgCor, y = -log10(q))) + geom_point(aes(color = highlight), alpha = 0.6, show.legend = F) + 
  geom_text(data = subset(TF_enrich_sub, highlight), aes(label = TF), nudge_x = -0.01, nudge_y = 0.6) + 
  scale_color_manual(values = c("grey50", "tomato")) + 
  xlab("Average TF-gene expression correlation") + ylab(expression(paste(-Log[10], " (FDR)"))) + 
  coord_cartesian(ylim = c(0, 10), clip = "off") + 
  theme(aspect.ratio = 1)

dev.off()
