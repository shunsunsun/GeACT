# cell-cell communication
setwd("~/lustre/06-Human_cell_atlas/pooled_data/01_stomach/")

library("reshape2")
library("glue")
library("ggplot2")
suppressMessages(library("cowplot"))
library("ggtext")
suppressMessages(library("igraph"))
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
ctMetaData <- subset(ctMetaData, stage == "14w")[-4]

# color
load("../../pooled_data/All/03-expression/merged/cellCluster/ct_color.RData")
cg_color <- c("tomato", "orange", "#FB9A99", "#00BF74", "#8258FA",
              "#D0A9F5", "purple2", "#D358F7", "#419FDE", "#F2E100",
              ct_color["Ovary"], ct_color["Testis"], "firebrick3", "grey80", 
              "#8258FA", "skyblue", "#B3B3B3")
names(cg_color) <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "B", 
                     "DC/Macrophage", "NKT", "T", "Glial", "FGC", "Granulosa", "Sertoli", "Erythrocyte", "Other", 
                     "Immune", "CACNA1A", "Unknown")

# load cellPhoneDB result
all_pval <- read.table(file.path(OUT, "pvalues_fixed.txt"), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
all_means <- read.table(file.path(OUT, "means_fixed.txt"), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
all(all_pval$id_cp_interaction == all_means$id_cp_interaction) & all(colnames(all_pval) == colnames(all_means))

# check receptor-ligand order
as.data.frame(table(all_pval[, c("receptor_a", "receptor_b")]))
# extract receptor-ligand pair
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
# write
write.table(x = all_pm, file = file.path(OUT, "cellComm_raw_St_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = subset(all_pm, pvalue < 1e-5 & mean >= 0.5)[, c(1, 5:16)], file = file.path(OUT, "cellComm_sig_St_byGroup.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
#

# 1. select major cell groups ----
ct_sel <- c("Epithelial", "Endothelial", "Fibroblast", "Glial", "Erythrocyte")
ctPair_sel <- apply(data.table::CJ(ct_sel, ct_sel, sorted = F), 1, function(x) { paste(x[1], x[2]) })
all_pm_sub <- subset(all_pm, ct_a %in% ct_sel & ct_b %in% ct_sel)
# filtering
all_pm_ftd <- subset(all_pm_sub, (is_integrin == "False") & (! is.na(gene_b_new)) & pvalue <= 1e-5 & mean >= 0.5)
# rescue terms
all_pm_ftr <- subset(all_pm_sub, genePair %in% all_pm_ftd$genePair & pvalue != 1)
all_pm_ftr$ctPair <- factor(all_pm_ftr$ctPair, levels = ctPair_sel)
all_pm_ftr$ctPair_fmt <- factor(all_pm_ftr$ctPair_fmt, levels = unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[match(levels(all_pm_ftr$ctPair), unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[, "ctPair"]), "ctPair_fmt"])

# heatmap
ga <- ggplot(all_pm_ftr, aes(x = ctPair_fmt, y = genePair_fmt)) + geom_point(aes(size=-log10(pvalue),color=mean)) + #shape = "\u25D6"
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") + 
  geom_vline(xintercept = 1:4 * 5 + 0.5, linetype = "dashed", color = "grey60") + 
  scale_color_gradientn(colors=colorRampPalette(c("grey80", "yellow", "red"), alpha=TRUE)(n=100)) + 
  labs(color = parse(text = "Log[2]~(CPM+1)"), size = parse(text = "-Log[10]~(p-value)")) + 
  theme(axis.text.x = element_markdown(size = 10, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_markdown(size = 10), 
        axis.title = element_blank(), 
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(margin = margin(r = 0.25, unit = "cm")), 
        legend.position = "bottom", legend.justification = c(0.5, 0.5), legend.margin = margin(t = -0.2, l = 0.2, r = 0.2, unit = "cm")) + 
  guides(color = guide_colorbar(title.position = "top", barwidth = unit(5, units = "cm"), title.hjust = 0.5, order = 1), 
         size = guide_legend(title.position = "top", keywidth = unit(0.1, units = "cm"), title.hjust = 0.5, order = 2, label.theme = element_text(margin = margin(l = -0.1, unit = "cm"))))
ggsave(file.path(OUT, "cellComm_heatmap.pdf"), plot = ga, width = 7, height = 35, device = cairo_pdf)

# 2. select major cell groups and a part of gene pairs ----
ct_sel <- c("Epithelial", "Endothelial", "Fibroblast", "Glial", "Erythrocyte")
ctPair_sel <- apply(data.table::CJ(ct_sel, ct_sel, sorted = F), 1, function(x) { paste(x[1], x[2]) })
all_pm_sub <- subset(all_pm, ct_a %in% ct_sel & ct_b %in% ct_sel)
# filtering
all_pm_ftd <- subset(all_pm_sub, (is_integrin == "False") & (! is.na(gene_b_new)) & pvalue <= 1e-5 & mean >= 0.5)
# rescue terms
all_pm_ftr <- subset(all_pm_sub, genePair %in% all_pm_ftd$genePair & pvalue != 1)
all_pm_ftr$ctPair <- factor(all_pm_ftr$ctPair, levels = ctPair_sel)
all_pm_ftr$ctPair_fmt <- factor(all_pm_ftr$ctPair_fmt, levels = unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[match(levels(all_pm_ftr$ctPair), unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[, "ctPair"]), "ctPair_fmt"])
# select specific gene pairs
all_pm_sel <- subset(all_pm_ftr, (grepl("BMP|CXCL|FAM3C|IGF|VEGFA|WNT5A", gene_a_new) | grepl("NOTCH", gene_b_new)) & ! genePair %in% c("NOV NOTCH1", "WNT4 NOTCH1") & gene_a_new != "DLK1")

gp <- ggplot(all_pm_sel, aes(x = ctPair_fmt, y = genePair_fmt)) + geom_point(aes(size=-log10(pvalue),color=mean)) + #shape = "\u25D6"
  scale_x_discrete(drop = F) + 
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") + 
  geom_vline(xintercept = 1:4 * 5 + 0.5, linetype = "dashed", color = "grey60") + 
  scale_color_gradientn(colors=colorRampPalette(c("grey80", "yellow", "red"), alpha=TRUE)(n=100)) + 
  labs(color = parse(text = "Log[2]~(CPM+1)"), size = parse(text = "-Log[10]~(p-value)")) + 
  theme(axis.text.x = element_markdown(size = 10, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_markdown(size = 10), 
        axis.title = element_blank(), 
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(margin = margin(r = 0.25, unit = "cm")), 
        legend.position = "bottom", legend.justification = c(0.5, 0.5)) + 
  guides(color = guide_colorbar(title.position = "top", barwidth = unit(5, units = "cm"), title.hjust = 0.5, order = 1), 
         size = guide_legend(title.position = "top", keywidth = unit(0.1, units = "cm"), title.hjust = 0.5, order = 2, label.theme = element_text(margin = margin(l = -0.1, unit = "cm"))))
ggsave(file.path(OUT, "cellComm_heatmap_sel.pdf"), plot = gp, width = 7, height = 12)

# 3. use Fibro as the key cell type ----
ct_sel <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "Glial", "B", "DC/Macrophage", "T", "CACNA1A", "Erythrocyte")
ctPair_sel <- apply(data.table::CJ(ct_sel, ct_sel, sorted = F), 1, function(x) { paste(x[1], x[2]) })
all_pm_sub <- subset(all_pm, ct_a %in% ct_sel & ct_b %in% ct_sel & (ct_a == "Fibroblast" | ct_b == "Fibroblast"))
# filtering
all_pm_ftd <- subset(all_pm_sub, (is_integrin == "False") & (! is.na(gene_b_new)) & pvalue <= 1e-5 & mean >= 0.5)
# rescue terms
all_pm_ftr <- subset(all_pm_sub, genePair %in% all_pm_ftd$genePair & pvalue != 1)
all_pm_ftr$ctPair <- factor(all_pm_ftr$ctPair, levels = c(setdiff(c(grep("^Fibroblast", ctPair_sel, value = T), grep("Fibroblast$", ctPair_sel, value = T)), "Fibroblast Fibroblast"), "Fibroblast Fibroblast"))
all_pm_ftr$ctPair_fmt <- factor(all_pm_ftr$ctPair_fmt, levels = unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[match(levels(all_pm_ftr$ctPair), unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[, "ctPair"]), "ctPair_fmt"])

# heatmap
gf <- ggplot(all_pm_ftr, aes(x = ctPair_fmt, y = genePair_fmt)) + geom_point(aes(size=-log10(pvalue),color=mean)) + #shape = "\u25D6"
  scale_x_discrete(drop = T) + 
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") + 
  geom_vline(xintercept = (length(ct_sel) - 1) * 1:2 + 0.5, linetype = "dashed", color = "grey60") + 
  scale_color_gradientn(colors=colorRampPalette(c("grey80", "yellow", "red"), alpha=TRUE)(n=100)) + 
  labs(color = parse(text = "Log[2]~(CPM+1)"), size = parse(text = "-Log[10]~(p-value)")) + 
  theme(axis.text.x = element_markdown(size = 10, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_markdown(size = 10), 
        axis.title = element_blank(), 
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(margin = margin(r = 0.25, unit = "cm")), 
        legend.position = "bottom", legend.justification = c(0.5, 0.5), legend.margin = margin(t = -0.2, l = 0.2, r = 0.2, unit = "cm")) + 
  guides(color = guide_colorbar(title.position = "top", barwidth = unit(5, units = "cm"), title.hjust = 0.5, order = 1), 
         size = guide_legend(title.position = "top", keywidth = unit(0.1, units = "cm"), title.hjust = 0.5, order = 2, label.theme = element_text(margin = margin(l = -0.1, unit = "cm"))))
ggsave(file.path(OUT, "cellComm_heatmap_Fibro.pdf"), plot = gf, width = 8, height = 35, device = cairo_pdf)

# 4. use Fibro as the key cell type (cell type specific interaction) ----
ct_sel <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "Glial", "B", "DC/Macrophage", "T", "CACNA1A", "Erythrocyte")
ctPair_sel <- apply(data.table::CJ(ct_sel, ct_sel, sorted = F), 1, function(x) { paste(x[1], x[2]) })
all_pm_sub <- subset(all_pm, ct_a %in% ct_sel & ct_b %in% ct_sel & (ct_a == "Fibroblast" | ct_b == "Fibroblast"))
# filtering
all_pm_ftd <- subset(all_pm_sub, (is_integrin == "False") & (! is.na(gene_b_new)) & pvalue <= 1e-5 & mean >= 0.5)
all_pm_ftd <- subset(all_pm_ftd, genePair %in% names(table(all_pm_ftd$genePair))[table(all_pm_ftd$genePair) == 1])
# only show part of gene pairs
all_pm_ftd_sub <- subset(all_pm_ftd, group == "Paracrine signaling")
all_pm_ftd_sub <- all_pm_ftd_sub[order(all_pm_ftd_sub$pvalue, - all_pm_ftd_sub$mean), ]
all_pm_ftd_sub <- do.call("rbind", lapply(split(all_pm_ftd_sub, all_pm_ftd_sub$ctPair), function(x) { head(x, 2) }))
rownames(all_pm_ftd_sub) <- NULL
all_pm_fsp <- rbind(subset(all_pm_ftd, group != "Paracrine signaling"), all_pm_ftd_sub)
all_pm_fsp$ctPair <- factor(all_pm_fsp$ctPair, levels = c(setdiff(c(grep("^Fibroblast", ctPair_sel, value = T), grep("Fibroblast$", ctPair_sel, value = T)), "Fibroblast Fibroblast"), "Fibroblast Fibroblast"))
all_pm_fsp <- all_pm_fsp[order(all_pm_fsp$ctPair), ]
# rescue terms
all_pm_ftr <- subset(all_pm_sub, genePair %in% all_pm_fsp$genePair & pvalue != 1)
all_pm_ftr$ctPair <- factor(all_pm_ftr$ctPair, levels = c(setdiff(c(grep("^Fibroblast", ctPair_sel, value = T), grep("Fibroblast$", ctPair_sel, value = T)), "Fibroblast Fibroblast"), "Fibroblast Fibroblast"))
all_pm_ftr$ctPair_fmt <- factor(all_pm_ftr$ctPair_fmt, levels = unique(all_pm_sub[, c("ctPair", "ctPair_fmt")])[match(levels(all_pm_ftr$ctPair), unique(all_pm_sub[, c("ctPair", "ctPair_fmt")])[, "ctPair"]), "ctPair_fmt"])
all_pm_ftr <- all_pm_ftr[order(all_pm_ftr$ctPair), ]
all_pm_ftr$genePair <- factor(all_pm_ftr$genePair, levels = rev(unique(all_pm_fsp$genePair)))
all_pm_ftr$genePair_fmt <- factor(all_pm_ftr$genePair_fmt, levels = unique(all_pm_ftr[, c("genePair", "genePair_fmt")])[match(levels(all_pm_ftr$genePair), unique(all_pm_ftr[, c("genePair", "genePair_fmt")])[, "genePair"]), "genePair_fmt"])

# heatmap
gs <- ggplot(all_pm_ftr, aes(x = ctPair_fmt, y = genePair_fmt)) + geom_point(aes(size=-log10(pvalue),color=mean)) + #shape = "\u25D6"
  scale_x_discrete(drop = F) + 
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") + 
  geom_vline(xintercept = (length(ct_sel) - 1) * 1:2 + 0.5, linetype = "dashed", color = "grey60") + 
  scale_color_gradientn(colors=colorRampPalette(c("grey80", "yellow", "red"), alpha=TRUE)(n=100)) + 
  labs(color = parse(text = "Log[2]~(CPM+1)"), size = parse(text = "-Log[10]~(p-value)")) + 
  theme(axis.text.x = element_markdown(size = 10, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_markdown(size = 10), 
        axis.title = element_blank(), 
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(margin = margin(r = 0.25, unit = "cm")), 
        legend.position = "bottom", legend.justification = c(0.5, 0.5), legend.margin = margin(t = -0.2, l = 0.2, r = 0.2, unit = "cm")) + 
  guides(color = guide_colorbar(title.position = "top", barwidth = unit(5, units = "cm"), title.hjust = 0.5, order = 1), 
         size = guide_legend(title.position = "top", keywidth = unit(0.1, units = "cm"), title.hjust = 0.5, order = 2, label.theme = element_text(margin = margin(l = -0.1, unit = "cm"))))
ggsave(file.path(OUT, "cellComm_heatmap_Fibro_spe.pdf"), plot = gs, width = 7, height = 10, device = cairo_pdf)

# 5. network plot (Fibroblast) ----
ct_sel <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "Glial", "B", "DC/Macrophage", "T", "CACNA1A", "Erythrocyte")
ctPair_sel <- apply(data.table::CJ(ct_sel, ct_sel, sorted = F), 1, function(x) { paste(x[1], x[2]) })
all_pm_sub <- subset(all_pm, ct_a %in% ct_sel & ct_b %in% ct_sel & (ct_a == "Fibroblast" | ct_b == "Fibroblast"))
# filtering
all_pm_ftd <- subset(all_pm_sub, (is_integrin == "False") & (! is.na(gene_b_new)) & pvalue <= 1e-5 & mean >= 0.5)
# no rescue terms
all_pm_ftr <- all_pm_ftd
all_pm_ftr$ctPair <- factor(all_pm_ftr$ctPair, levels = c(setdiff(c(grep("^Fibroblast", ctPair_sel, value = T), grep("Fibroblast$", ctPair_sel, value = T)), "Fibroblast Fibroblast"), "Fibroblast Fibroblast"))
all_pm_ftr$ctPair_fmt <- factor(all_pm_ftr$ctPair_fmt, levels = unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[match(levels(all_pm_ftr$ctPair), unique(all_pm_ftr[, c("ctPair", "ctPair_fmt")])[, "ctPair"]), "ctPair_fmt"])
# gene pair as labels
all_pm_fsm <- do.call("rbind", lapply(split(all_pm_fsp, all_pm_fsp$ctPair), function(x) { head(x, 1) }))
rownames(all_pm_fsm) <- NULL
subset(all_pm_fsm, ct_a == "Fibroblast")

cellComm <- data.frame(table(all_pm_ftr$ctPair), stringsAsFactors = F)
colnames(cellComm) <- c("ctPair", "nGenePair")
cellComm <- subset(cellComm, nGenePair > 0)
cellComm <- merge(cellComm, unique(all_pm[, c("ctPair", "ct_a", "ct_b")]), by = "ctPair", sort = F)
cellComm <- merge(cellComm, subset(all_pm_fsm, ct_a == "Fibroblast" | ct_b == "Fibroblast", c("ctPair", "gene_a_new", "gene_b_new", "genePair")), by = "ctPair", all.x = T, sort = F)
colnames(cellComm)[5] <- "label"
nodes <- data.frame(ident = ct_sel, stringsAsFactors = F)
edges <- cellComm[, c(3,4,2,1,5)]

# node attr
nodes$color <- cg_color[nodes$ident]
# nodes$label <- nodes$gene
# nodes$label.cex <- 0.4
nodes$frame.color <- "grey90"
# nodes$frame.width <- 1

# edge attr
ecols <- c("sandybrown", "palevioletred1", "#00BF74")
edges$width <- 1
edges$lty <- "solid"
edges$color <- "grey90"
edges$color[edges$ct_a == "Fibroblast" & ! is.na(edges$label)] <- ecols[1]
edges$color[edges$ct_b == "Fibroblast" & ! is.na(edges$label)] <- ecols[2]
edges$color[edges$ct_a == "Fibroblast" & edges$ct_b == "Fibroblast"] <- ecols[3]

g <- graph_from_data_frame(edges, directed = T, vertices = nodes)
if(grepl("_14w", strsplit(getwd(), split = "/")[[1]][6])) {
  print("Use the node positions in 20w.")
  vertex <- read.table(file.path("../../pooled_data", ts_ip, OUT, "vertex.txt"), header = F, sep = "\t", stringsAsFactors = F)[, 1]
  coords <- as.matrix(read.table(file.path("../../pooled_data", ts_ip, OUT, "coords.txt"), header = F, sep = "\t", stringsAsFactors = F))
  coords <- coords[match(as_ids(V(g)), vertex), ]
} else {
  print("Calculate node positions.")
  set.seed(1234)
  #coords <- layout_with_fr(g, niter = 10000)
  coords <- layout_as_star(g, center = "Fibroblast")
  rownames(coords) <- as_ids(V(g))
}

edge.label.attr <- eLableAttr(g, coords, ml_ratio = 0.4, ol_ratio = 0.1, loop_adj_x = 0.8, loop_adj_y = 0.25)
edge.label.attr["CACNA1A|Fibroblast", 1] <- edge.label.attr["CACNA1A|Fibroblast", 1] - 0.03
edge.label.attr["Fibroblast|Glial", 1] <- edge.label.attr["Fibroblast|Glial", 1] - 0.05
edge.label.attr["Fibroblast|T", 1] <- edge.label.attr["Fibroblast|T", 1] - 0.03
edge.label.attr["T|Fibroblast", 1] <- edge.label.attr["T|Fibroblast", 1] - 0.03
edge.label.attr$color <- ifelse(grepl("^Fibroblast", rownames(edge.label.attr)), ecols[1], ecols[2])
edge.label.attr["Fibroblast|Fibroblast", "color"] <- ecols[3]

pdf(file.path(OUT, "cellComm_network_Fibro_one_other.pdf"), width = 5, height = 5)

par(mar = c(0,0,0,0))
plot.igraph_new(g, layout = coords, asp = 1, 
                vertex.size = 12, #vertex.size = 6, sqrt(degree(g)) * 3
                vertex.shape = "circle", 
                #vertex.frame.color = "white", 
                vertex.label.color = "black", 
                #vertex.label.cex = 0.2, 
                edge.curved = 0.1, 
                edge.arrow.size= 0.5, 
                edge.label.x = edge.label.attr$x, 
                edge.label.y = edge.label.attr$y, 
                #edge.label.cex = 0.8, 
                edge.label.color = edge.label.attr$color, 
                edge.label.srt = edge.label.attr$srt, 
                xlim = c(-1.2, 1.1) #, ylim = c(-0.9, 0.9)
)
par(mar = c(5,4,4,2) + 0.1)

dev.off()

# write structure
write.table(x = nodes, file = file.path(OUT, "nodes.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = edges, file = file.path(OUT, "edges.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = as_ids(V(g)), file = file.path(OUT, "vertex.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(x = coords, file = file.path(OUT, "coords.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# X. save info ----
save.image(file = paste0(OUT, "/cellComm.RData"))
