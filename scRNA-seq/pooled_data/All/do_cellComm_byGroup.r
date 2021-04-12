# cell-cell communication
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("reshape2")
library("glue")
library("ggplot2")
suppressMessages(library("cowplot"))
#library("ggtext")
#suppressMessages(library("igraph"))
suppressMessages(library("ComplexHeatmap"))
source("../../scripts/cellComm_tools.r")

samplingPos <- "byGroup"
OUT <- paste0("03-expression/merged/cellComm/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/cellComm.RData"))

# load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Genomes/human/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(geneID) <- c("ensembl", "symbol")
geneID$ensembl_alt <- gsub("\\.[0-9]+", "", geneID$ensembl)

# load cell type meta
ts_ip <- gsub(".*/", "", getwd())
ctMetaData <- read.table(file = file.path("../../pooled_data_all/", ts_ip, "cellType_metatable.txt"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")

# load cellPhoneDB result
cellComm_LS <- lapply(list.files("..", pattern = "^[0-9]"), function(x) {
  cat(">", x, "\n")
  all_pval <- read.table(file = paste0("../", x, "/03-expression/merged/cellComm/byGroup/pvalues_fixed.txt"), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  all_means <- read.table(file = paste0("../", x, "/03-expression/merged/cellComm/byGroup/means_fixed.txt"), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  all(all_pval$id_cp_interaction == all_means$id_cp_interaction) & all(colnames(all_pval) == colnames(all_means))
  
  # check ligand-receptor order
  as.data.frame(table(all_pval[, c("receptor_a", "receptor_b")]))
  # extract ligand-receptor pair
  all_pval <- subset(all_pval, receptor_a == "False" & receptor_b == "True")
  all_means <- subset(all_means, receptor_a == "False" & receptor_b == "True")
  # remove duplicated terms
  all_pval <- all_pval[! duplicated(all_pval$interacting_pair), ]
  all_means <- all_means[! duplicated(all_means$interacting_pair), ]
  # combine
  all_pval_melted <- melt(all_pval, id.vars = colnames(all_pval)[1:11], variable.name = "ctPair", value.name = "pvalue")
  all_means_melted <- melt(all_means, id.vars = colnames(all_means)[1:11], variable.name = "ctPair", value.name = "mean")
  all_pm <- data.frame(all_pval_melted, mean = log2(all_means_melted$mean + 1), stringsAsFactors = F)
  # gene ID mapping
  all_pm$gene_a_new <- geneID$symbol[match(all_pm$gene_a, geneID$ensembl_alt)]
  all_pm$gene_b_new <- geneID$symbol[match(all_pm$gene_b, geneID$ensembl_alt)]
  all_pm$genePair <- paste(all_pm$gene_a_new, all_pm$gene_b_new, sep = " ")
  # parse cellType Pair
  all_pm$ct_a <- gsub("\\|.*", "", all_pm$ctPair)
  all_pm$ct_b <- gsub(".*\\|", "", all_pm$ctPair)
  all_pm$ctPair <- gsub("|", " ", all_pm$ctPair, fixed = T)
  # set color
  all_pm$genePair_fmt <- glue("<span style='color:#0072B2'>{all_pm$gene_a_new}</span> <span style='color:#009E73'>{all_pm$gene_b_new}</span>")
  all_pm$ctPair_fmt <- glue("<span style='color:#0072B2'>{all_pm$ct_a}</span> <span style='color:#009E73'>{all_pm$ct_b}</span>")
  # set group by if secreted
  all_pm$group <- ifelse(all_pm$secreted == "True", "Paracrine signaling", "Juxtacrine signaling")
  # tune values for plot
  all_pm$pvalue[all_pm$pvalue == 0] <- 1e-10
  # add tissue info
  all_pm$tissue <- gsub("_", " ", gsub("^[0-9]+_", "", x))
  return(all_pm)
})
all_pm <- do.call("rbind", cellComm_LS)
# sig only
all_pm_sig <- subset(all_pm, pvalue < 1e-5 & mean >= 0.5)
length(unique(paste(all_pm_sig$tissue, all_pm_sig$id_cp_interaction, sep = ".")))
length(unique(all_pm_sig$id_cp_interaction))
# write
write.table(x = all_pm, file = file.path(OUT, "cellComm_raw_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = all_pm_sig[, c(1, 5:11, 13:16, 18:19, 23)], file = file.path(OUT, "cellComm_sig_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#

# 1. select major cell groups ----
ct_sel <- c("Epithelial", "Endothelial", "Fibroblast", "Glial", "Erythrocyte")
ts_sel <- intersect(unique(ctMetaData$tissue), c("stomach", "small intestine", "kidney", "lung", "pancreas", "spleen", "testis", "liver", "ovary"))
ctPair_sel <- apply(data.table::CJ(ct_sel, ct_sel, sorted = F), 1, function(x) { paste(x[1], x[2]) })
all_pm_sub <- subset(all_pm, ct_a %in% ct_sel & ct_b %in% ct_sel & tissue %in% ts_sel)
# filtering
all_pm_ftd <- subset(all_pm_sub, (is_integrin == "False") & (! is.na(gene_b_new)) & pvalue <= 1e-5 & mean >= 0.5)
# rescue terms
#all_pm_ftr <- subset(all_pm_sub, genePair %in% all_pm_ftd$genePair & pvalue != 1)
all_pm_ftr <- all_pm_ftd
all_pm_ftr$ctPair <- factor(all_pm_ftr$ctPair, levels = ctPair_sel)
all_pm_ftr$ctPair_fmt <- factor(all_pm_ftr$ctPair_fmt, levels = unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[match(levels(all_pm_ftr$ctPair), unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[, "ctPair"]), "ctPair_fmt"])

all_pm_acs <- acast(data = all_pm_ftr, formula = genePair ~ tissue + ctPair, value.var = "pvalue")
all_pm_acs <- -log10(all_pm_acs)
all_pm_acs[! is.na(all_pm_acs)] <- 1
all_pm_acs[is.na(all_pm_acs)] <- 0
annotation_col <- data.frame(Tissue = gsub("_.*", "", colnames(all_pm_acs)), ctPair = gsub(".*_", "", colnames(all_pm_acs)))
annotation_col$sendingCell <- gsub(" .*", "", annotation_col$ctPair)
annotation_col$receivingCell <- gsub(".* ", "", annotation_col$ctPair)
rownames(annotation_col) <- colnames(all_pm_acs)
#pheatmap::pheatmap(all_pm_acs, annotation_col = annotation_col[, -2])

load("03-expression/merged/cellCluster/ct_color.RData")
ccmap_ts <- Hmisc::capitalize(gsub("_.*", "", colnames(all_pm_acs)))
ccmap_ts_color <- ct_color[ccmap_ts]
names(ccmap_ts_color) <- ccmap_ts
ccmap_ts <- factor(ccmap_ts, levels = Hmisc::capitalize(ts_sel))

load("03-expression/merged/cellCluster/cg_color.RData")
ccmap_cta <- gsub(" .*", "", gsub(".*_", "", colnames(all_pm_acs)))
ccmap_cta_color <- cg_color[ccmap_cta]
ccmap_cta <- factor(ccmap_cta, levels = ct_sel)
ccmap_ctb <- gsub(".* ", "", gsub(".*_", "", colnames(all_pm_acs)))
ccmap_ctb_color <- cg_color[ccmap_ctb]
ccmap_ctb <- factor(ccmap_ctb, levels = ct_sel)

pdf(paste0(OUT, "/cellComm_map.pdf"), width = 8.25, height = 4.5)

ht <- Heatmap(matrix = all_pm_acs, name = "Communication", col = colorRampPalette(c("white", "limegreen"))(100),
              border = "grey50", 
              cluster_rows = T, cluster_columns = T, show_parent_dend_line = F, show_row_names = F, show_column_names = F, 
              #row_split = factor(cellType_sub$group, levels = unique(cellType_sub$group)), column_split = riskGene_sub$diseaseName_new, 
              row_title_gp = gpar(size = 12), row_title_rot = 0, 
              top_annotation = columnAnnotation(Organ = ccmap_ts, Sending = ccmap_cta, Receiving = ccmap_ctb, 
                                                col = list(Organ = ccmap_ts_color, Sending = ccmap_cta_color, Receiving = ccmap_ctb_color), 
                                                annotation_name_gp = gpar(fontsize = 12), show_legend = F), 
              show_heatmap_legend = F)

lgd1 <- Legend(at = Hmisc::capitalize(ts_sel), title = "Organ", legend_gp = gpar(fill = ct_color[Hmisc::capitalize(ts_sel)]), title_gp = gpar(fontsize = 12), title_gap = unit(2, "mm"), labels_gp = gpar(fontsize = 12))
lgd2 <- Legend(at = ct_sel, title = "Cell group", legend_gp = gpar(fill = cg_color[ct_sel]), title_gp = gpar(fontsize = 12), title_gap = unit(2, "mm"), labels_gp = gpar(fontsize = 12))
lgd <- packLegend(lgd1, lgd2, row_gap = unit(0.6, "cm"))
draw(ht, row_title = "Cellular communications                    ", padding = unit(c(5.5, 5.5, 5.5, 60), units = "points"))
draw(lgd, x = unit(0.9, "npc"), y = unit(0.375, "npc"), just = c("center", "center"))

dev.off()

# stat
nCtPair <- rowSums(all_pm_acs)
cellComm_byOrgan <- sapply(ts_sel, function(x) {
  all_pm_acs_sel <- all_pm_acs[, grep(paste0("^", x), colnames(all_pm_acs)), drop = F]
  y <- rowSums(all_pm_acs_sel)
  return(y)
})
nOrgan <- apply(cellComm_byOrgan > 0, 1, sum)
all(names(nCtPair) == names(nOrgan))
cellCommStat <- data.frame(nCtPair, nOrgan, stringsAsFactors = F)
cellCommStat$genePair <- rownames(cellCommStat)

# ovary specific
cellComm_ovary_spec <- cellComm_byOrgan[rowSums(cellComm_byOrgan[, c(1:7, 9)]) == 0 & cellComm_byOrgan[, 8] > 0, , drop = F]
cellComm_ovary_spec <- cellComm_ovary_spec[order(cellComm_ovary_spec[, 8], decreasing = T), ]
cellComm_ovary_spec
# testis specific
cellComm_testis_spec <- cellComm_byOrgan[rowSums(cellComm_byOrgan[, c(1:7, 8)]) == 0 & cellComm_byOrgan[, 9] > 0, , drop = F]
cellComm_testis_spec <- cellComm_testis_spec[order(cellComm_testis_spec[, 8], decreasing = T), ]
cellComm_testis_spec

cellCommStat$color <- "grey"
cellCommStat$color[cellCommStat$nOrgan == 9] <- "cornflowerblue"
cellCommStat$color[cellCommStat$genePair %in% rownames(cellComm_ovary_spec)[1]] <- ct_color["Ovary"]
cellCommStat$color[cellCommStat$genePair %in% rownames(cellComm_testis_spec)[1]] <- ct_color["Testis"]
cellCommStat <- cellCommStat[c(which(cellCommStat$color == "grey"), which(cellCommStat$color != "grey")), ]

pdf(paste0(OUT, "/cellComm_stat.pdf"), width = 4.5, height = 4.5)

set.seed(1)
ggplot(cellCommStat, aes(x = nOrgan, y = nCtPair)) + geom_jitter(aes(color = I(color)), width = 0.1) + 
  #geom_bin2d() + 
  ggrepel::geom_text_repel(data = subset(cellCommStat, nOrgan == 9), aes(label = genePair), color = "cornflowerblue", point.padding = 0.5, nudge_y = 7.5) + 
  ggrepel::geom_text_repel(data = cellCommStat[c(rownames(cellComm_ovary_spec)[1], rownames(cellComm_testis_spec)[1]), ], aes(label = genePair), 
                           color = ct_color[c("Ovary", "Testis")], point.padding = 0.5, nudge_x = 1, nudge_y = 12) + 
  xlab("Organ number") + ylab("Cell type pair number") + 
  scale_x_continuous(limits = c(0, 9.08), breaks = seq(0, 9, by = 3)) + 
  theme(aspect.ratio = 1)

dev.off()

# X. save info ----
save.image(file = paste0(OUT, "/cellComm.RData"))
