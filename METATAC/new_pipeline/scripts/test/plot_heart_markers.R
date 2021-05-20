setwd("/data/Lab/otherwork/GeACT/ATAC/data/19-22w/13_heart/results/")

marker_list <- read.table("/data/Lab/otherwork/GeACT/ATAC/meta/1-s2.0-S2211124719301081-mmc3.txt", 
           header = T, sep = "\t", stringsAsFactors = F)

proj <- loadArchRProject("ArchR/ArchR_output/", showLogo = F)
geneScore <- getMatrixFromProject(proj)
geneScore <- summarizedExperiment2Seurat(se = geneScore, assay = "ACTIVITY")

top5 <- marker_list %>% group_by(Cluster) %>% top_n(n = 5, wt = myAUC)

pdf("19_22w_heart_markers.pdf")
plotMarker(proj, geneScore, genes = top5$Gene)
dev.off()

# "EPCAM", "KRT19", "PTPRC", 
pdf("test.pdf")
plotMarker(proj, geneScore, genes = c("CD19", "CD79A", "MACRO", "HBG1", "HBG2", "HBB"))
dev.off()

pdf("test.pdf")
plotMarker(proj, geneScore, genes = c("CD3D", "CD3E", "CD4", "CD8A"))
dev.off()

pdf("test.pdf")
plotMarker(proj, geneScore, genes = c("PLP1", "S100P", "MARCO"))
dev.off()

pdf("test.pdf")
plotMarker(proj, geneScore, genes = c("ACTG2", "MYB"))
dev.off()
