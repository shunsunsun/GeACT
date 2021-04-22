# peak to gene link
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/02_small_intestine/")

suppressMessages({
  library("Seurat", lib.loc = "/rd1/apps/R-3.5.1/library_ext/")
  library("ggplot2")
  library("cowplot")
})
source("../../scripts/peak2gene_tools.r")

OUT <- "05-integration/merged/peak2gene"
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/peak2gene.RData"))

ts_ip <- gsub(".*/", "", getwd())

# load gene pos
gene_pos <- read.table(file = "../../Genomes/human/gene_pos.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_pos) <- c("ensembl", "gene", "chr", "left", "right", "strand", "tss")

# load links
link_11_14 <- read.table(file = file.path("../../ATAC/data/paired", ts_ip, "peak2gene/11-14w/Peak2GeneLinks/peak2gene.txt"), header = T, sep = "\t", stringsAsFactors = F)
link_19_22 <- read.table(file = file.path("../../ATAC/data/paired", ts_ip, "peak2gene/19-22w/Peak2GeneLinks/peak2gene.txt"), header = T, sep = "\t", stringsAsFactors = F)
link_11_14_ftd <- subset(link_11_14, abs(Correlation) >= 0.3)
link_19_22_ftd <- subset(link_19_22, abs(Correlation) >= 0.3)
nPeak_ftd <- merge(as.data.frame(table(link_11_14_ftd$gene)), as.data.frame(table(link_19_22_ftd$gene)), by = "Var1", sort = F)
colnames(nPeak_ftd) <- c("gene", "nPeak_ftd.x", "nPeak_ftd.y")
nPeak_ftd$gene <- as.character(nPeak_ftd$gene)

# load peaks
# peak_11_14 <- readRDS(file = file.path("../../ATAC/data/paired", ts_ip, "peak2gene/11-14w/Peak2GeneLinks/seuratATAC-Group-KNN.rds"))
# peak_19_22 <- readRDS(file = file.path("../../ATAC/data/paired", ts_ip, "peak2gene/19-22w/Peak2GeneLinks/seuratATAC-Group-KNN.rds"))

# differentially opened peaks
# peak_merged <- merge(x = peak_11_14, y = peak_19_22, add.cell.ids = c("11-14w", "19-22w"), project = "merged")
# peak_merged@meta.data$stage <- gsub("w_.*", "w", rownames(peak_merged@meta.data))
# #VlnPlot(peak_merged, features = c("nFeature_ATAC", "nCount_ATAC"), group.by = "stage")
# peak_merged <- NormalizeData(peak_merged, normalization.method = "LogNormalize", scale.factor = 10000)
# peak_merged <- ScaleData(peak_merged, features = rownames(peak_merged))
# markPeaks_11_14 <- FindMarkers(peak_merged, ident.1 = "11-14w", ident.2 = "19-22w", min.pct = 0.25, only.pos = T, group.by = "stage")
# markPeaks_19_22 <- FindMarkers(peak_merged, ident.1 = "19-22w", ident.2 = "11-14w", min.pct = 0.25, only.pos = T, group.by = "stage")
# markPeaks_11_14_ftd <- subset(markPeaks_11_14, p_val_adj < 0.05)
# markPeaks_19_22_ftd <- subset(markPeaks_19_22, p_val_adj < 0.05)

# load expr cmp results
exprCmp <- read.table(file = "03-expression/merged/cellCluster/Seurat_expr.markers_byGroup.txt", header = T, sep = "\t", stringsAsFactors = F)
exprCmp <- subset(exprCmp, cluster == "Epithelial")
exprCmp_nonDE <- subset(exprCmp, p_val > 0.75 & abs(avg_logFC) < log(1.05))
dim(exprCmp_nonDE)

# differential sig links
link <- merge(link_11_14, link_19_22, by = c("peak", "gene"), sort = F)
link$diff_abs <- abs(link$Correlation.y) - abs(link$Correlation.x)
dsl_a <- subset(link, diff_abs < -0.3 & FDR.x < 1e-5 & FDR.y > 0.05)
dsl_b <- subset(link, diff_abs > 0.3 & FDR.x > 0.05 & FDR.y < 1e-5)
dsl_s <- rbind(dsl_a, dsl_b)
dsl_s$peakChr <- gsub(":.*", "", dsl_s$peak)
dsl_s$peakLeft <- as.numeric(gsub("-.*", "", gsub(".*:", "", dsl_s$peak)))
dsl_s$peakRight <- as.numeric(gsub(".*-", "", gsub(".*:", "", dsl_s$peak)))
dsl_s$peakCenter <- (dsl_s$peakLeft + dsl_s$peakRight) / 2
nPeak_dif <- merge(as.data.frame(table(dsl_a$gene)), as.data.frame(table(dsl_b$gene)), by = "Var1", sort = F)
colnames(nPeak_dif) <- c("gene", "nPeak_dif.x", "nPeak_dif.y")
nPeak_dif$gene <- as.character(nPeak_dif$gene)
nPeak_dif <- merge(nPeak_dif, nPeak_ftd, by = "gene", sort = F)
nPeak_dif$nonDE <- nPeak_dif$gene %in% exprCmp_nonDE$gene
nPeak_dif_nonDE <- subset(nPeak_dif, nonDE)
nPeak_dif_nonDE_11 <- subset(nPeak_dif_nonDE, nPeak_dif.x == 1 & nPeak_dif.y == 1)
cat("Could be used:", nrow(nPeak_dif_nonDE_11), "/", nrow(nPeak_dif_nonDE), "/", nrow(nPeak_dif), "\n")

# plot
pdf(file = paste0(OUT, "/dsl_ATF2.pdf"), width = 6, height = 5)
gp_LS <- do_plotTrack(case_gene = "ATF2", margin.right = 35)
plot_grid(gp_LS[[1]], NULL, gp_LS[[2]], NULL, gp_LS[[3]], ncol = 1, align = "hv", rel_heights = c(1.05, -0.2, 1, -0.2, 1))
dev.off()

pdf(file = paste0(OUT, "/dsl_TLK2.pdf"), width = 6, height = 5)
gp_LS <- do_plotTrack(case_gene = "TLK2", margin.right = 15)
plot_grid(gp_LS[[1]], NULL, gp_LS[[2]], NULL, gp_LS[[3]], ncol = 1, align = "hv", rel_heights = c(1.05, -0.2, 1, -0.2, 1))
dev.off()

pdf(file = paste0(OUT, "/dsl_nonDE_case.pdf"), width = 6, height = 3)

for(case_gene in nPeak_dif_nonDE_11$gene) {
  do_plotTrack(case_gene)
}

dev.off()

# X. save ----
#save.image(file = paste0(OUT, "/peak2gene.RData"))
