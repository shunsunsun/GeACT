Args <- commandArgs()
setwd(Args[6])
dir.create("Seurat_integration",recursive=T)

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(gridExtra)

# # Gene activity quantification
# peaks <- t(read.table(gzfile("insertcount_All_NarrowPeaks_split_filter.tsv.gz"),
#                     sep='\t',header=T,row.names = 1,stringsAsFactors = F, check.names = F)[,-1])
# # create a gene activity matrix from cicero
# activity.matrix <- read.table(gzfile("cicero/cicero_gene_activities_unnorm.tsv.gz"),
#                               sep='\t',header=T,row.names = 1,stringsAsFactors = F, check.names = F)
# 
# # Object setup
# expr.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "ATAC")
# expr.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
# # meta <- read.table("meta_info.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = F)
# # names(meta)
# # expr.atac <- AddMetaData(expr.atac, metadata = meta$ident,"Ident_old")
# # expr.atac <- subset(expr.atac, subset = nCount_ATAC > 5000)
# expr.atac$tech <- "atac"
# 
# pdf("Seurat_integration/Seurat_integration.pdf", width = 6, height = 6, useDingbats = F)
# 
# # Data preprocessing
# DefaultAssay(expr.atac) <- "ACTIVITY"
# # expr.atac <- FindVariableFeatures(expr.atac)
# expr.atac <- NormalizeData(expr.atac)
# expr.atac <- ScaleData(expr.atac,features = rownames(expr.atac[['ACTIVITY']]))
# print(dim(expr.atac@assays$ACTIVITY@scale.data))
# 
# DefaultAssay(expr.atac) <- "ATAC"
# VariableFeatures(expr.atac) <- names(which(Matrix::rowSums(expr.atac) > 2*ncol(expr.atac)))
# expr.atac <- RunLSI(expr.atac, n = 50, scale.max = NULL)
# expr.atac <- RunTSNE(expr.atac, reduction = "lsi",dims = 1:50)
# expr.atac <- NormalizeData(expr.atac)
# expr.atac <- ScaleData(expr.atac,features = rownames(expr.atac[['ATAC']]))
# print(dim(expr.atac@assays$ATAC@scale.data))
# 
# # Idents(expr.atac) <- expr.atac$Ident_old
# # print(DimPlot(expr.atac, group.by = "Ident_old", label = T, repel = T) +
# #           ggtitle("scATAC-seq cells") + NoLegend() + scale_colour_hue(drop = F))
# 
expr.rna <- readRDS(Args[7])
expr.rna[['ident']] <- Idents(expr.rna)
expr.rna$tech <- "rna"
# 
# transfer.anchors <- FindTransferAnchors(
#     reference = expr.rna, query = expr.atac,
#     features = VariableFeatures(object = expr.rna),
#     reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
# 
# transfer.anchors.table <- data.frame(transfer.anchors@anchors)
# transfer.anchors.table[,1] <- gsub("_reference$","",
#                                    transfer.anchors@reference.cells[transfer.anchors.table[,1]])
# transfer.anchors.table[,2] <- gsub("_query$","",
#                                    transfer.anchors@query.cells[transfer.anchors.table[,2]])
# names(transfer.anchors.table)[c(1,2)] <- c("RNA_cell","ATAC_cell")
# write.table(transfer.anchors.table,"Seurat_integration/transfer_anchors.txt",
#             sep='\t',quote = F,row.names = F,col.names = T)
# 
# celltype.predictions <- TransferData(
#     anchorset = transfer.anchors, refdata = expr.rna$ident,
#     weight.reduction = expr.atac[["lsi"]])
# expr.atac <- AddMetaData(expr.atac, metadata = celltype.predictions$predicted.id,"predicted.id")
# expr.atac <- AddMetaData(expr.atac, metadata = celltype.predictions$prediction.score.max,"prediction.score.max")
# # to make the colors match
# expr.atac$predicted.id <- factor(expr.atac$predicted.id, levels = levels(expr.rna))
# Idents(expr.atac) <- celltype.predictions$predicted.id
# write.table(celltype.predictions,"Seurat_integration/celltype_predictions.txt",
#             sep='\t',quote = F,row.names = T,col.names = T)
# 
# # Co-embedding
# # note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# # full transcriptome if we wanted to
# genes.use <- VariableFeatures(expr.rna)
# refdata <- GetAssayData(expr.rna, assay = "RNA", slot = "data")[genes.use, ]
# 
# # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
#                            weight.reduction = expr.atac[["lsi"]])
# 
# # this line adds the imputed data matrix to the expr.atac object
# expr.atac[["RNA"]] <- imputation
# saveRDS(expr.atac, file = "Seurat_integration/Seurat_expr_atac_integration.rds")
# coembed <- merge(x = expr.rna, y = expr.atac)
# 
# # Finally, we run PCA and TSNE on this combined object, to visualize the co-embedding of both
# # datasets
# coembed <- ScaleData(coembed, features = genes.use, do.scale = F)
# coembed <- RunPCA(coembed, features = genes.use, verbose = F)
# coembed <- RunTSNE(coembed, dims = 1:30)
# coembed$celltype <- ifelse(!is.na(coembed$ident), coembed$ident, coembed$predicted.id)
# saveRDS(coembed, file = "Seurat_integration/Seurat_expr_coembed.rds")

expr.atac <- readRDS("Seurat_integration/Seurat_expr_atac_integration.rds")
coembed <- readRDS("Seurat_integration/Seurat_expr_coembed.rds")

markers = c("EPCAM" = "EPCAM (CD326)","VIM" = "VIM","PTPRC" = "PTPRC (CD45)","HBG1" = "HBG1")
cellMetaData <- merge(expr.atac@meta.data, expr.atac@reductions$tsne@cell.embeddings, by = 0, sort = F)
cellMetaData <- merge(cellMetaData, t(expr.atac@assays$ACTIVITY@scale.data[names(markers),]), 
                      by.x = 'Row.names',by.y = 0, sort = F)
rownames(cellMetaData) <- cellMetaData$Row.names
cellMetaData <- cellMetaData[,-1]
write.table(cellMetaData,"Seurat_integration/meta_info_integration.txt",
            sep='\t',quote = F,row.names = T,col.names = T)

print(ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = orig.ident)) + geom_point(size = 0.5) +
          theme_cowplot() +
          theme(legend.justification = c("center")) +
          theme(aspect.ratio = 1, legend.box.margin = margin(t = -15), plot.margin = margin(t = 10, l = -20, r = 20)) +
          theme(legend.key.width = unit(0.2, "cm")) +
          labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") +
          guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))))

print(ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = predicted.id)) + geom_point(size = 0.5) +
          theme_cowplot() +
          theme(legend.justification = c("center")) +
          theme(aspect.ratio = 1, legend.box.margin = margin(t = -15), plot.margin = margin(t = 10, l = -20, r = 20)) +
          theme(legend.key.width = unit(0.2, "cm")) +
          labs(color = NULL) + xlab("t-SNE1") + ylab("t-SNE2") +
          guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.6))))

print(FeaturePlot(object = expr.atac, features = "prediction.score.max", reduction = "tsne",
                  pt.size = 0.5,label = T,repel = T))
print(DimPlot(expr.rna, group.by = "ident", label = T, repel = T) +
          ggtitle("scRNA-seq cells") + NoLegend())

print(DimPlot(coembed, group.by = "tech"))
print(DimPlot(coembed, group.by = "celltype", label = T, repel = T) + NoLegend())


DefaultAssay(expr.atac) <- "ACTIVITY"

plot_markers = list()

plot_markers[['EPCAM']] = ggplot(cellMetaData,aes(x=tSNE_1,y=tSNE_2)) +
    geom_point(aes(color=EPCAM),size=0.5,show.legend = F) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() + theme(aspect.ratio = 1) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "EPCAM (CD326)",x="",y="")

plot_markers[['VIM']] = ggplot(cellMetaData,aes(x=tSNE_1,y=tSNE_2)) +
    geom_point(aes(color=VIM),size=0.5,show.legend = F) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() + theme(aspect.ratio = 1) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "VIM",x="",y="")

plot_markers[['PTPRC']] = ggplot(cellMetaData,aes(x=tSNE_1,y=tSNE_2)) +
    geom_point(aes(color=PTPRC),size=0.5,show.legend = F) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() + theme(aspect.ratio = 1) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "PTPRC (CD45)",x="",y="")

plot_markers[['HBG1']] = ggplot(cellMetaData,aes(x=tSNE_1,y=tSNE_2)) +
    geom_point(aes(color=HBG1),size=0.5,show.legend = F) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() + theme(aspect.ratio = 1) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "HBG1",x="",y="")

print(grid.arrange(grobs=plot_markers,ncol=2, nrow=2, widths=c(3,3),heights=c(3,3)))

dev.off()

# expr.atac.filtered <- subset(expr.atac, subset = prediction.score.max > 0.5)
# 
# # Finding differentially expressed features (cluster biomarkers)
# pdf("Seurat_integration/Seurat_TSNE.pdf",width = 6, height = 6, useDingbats = F)
# 
# DefaultAssay(expr.atac.filtered) <- "ATAC"
# DimPlot(expr.atac.filtered, reduction = "tsne", group.by = "predicted.id", label = T, repel = T) +
#     NoLegend() + scale_colour_hue(drop = F)
# 
# markers <- FindAllMarkers(expr.atac.filtered, only.pos = T, min.pct = 0.1, logfc.threshold = 0.05)
# table(markers$cluster)
# # hist(-log10(markers$p_val_adj))
# # hist(markers$avg_logFC)
# markers_ftd <- markers[markers$p_val_adj <= 0.05,]
# table(markers_ftd$cluster)
# 
# top1 <- markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
# top5 <- markers_ftd %>% group_by(cluster) %>% top_n(5, avg_logFC)
# 
# VlnPlot(expr.atac.filtered, features = c("nCount_ATAC"), slot = "counts", log = T)
# VlnPlot(expr.atac.filtered, features = c("nFeature_ATAC"), slot = "counts", log = T)
# # VlnPlot(expr.atac.filtered, features = top1$gene, slot = "counts", log = T)
# for(i in top1$gene) {
#     print(FeaturePlot(object = expr.atac.filtered, features = i, reduction = "tsne",
#                       pt.size = 0.5,label = T,repel = T))
# }
# 
# write.table(markers_ftd,"Seurat_integration/cluster_marker_peak.txt",
#             row.names = T,col.names = T,sep='\t',quote=F)
# 
# DoHeatmap(expr.atac.filtered, features = top5$gene,disp.min = -3,disp.max = 3,size = 2.5)
# 
# 
# 
# DefaultAssay(expr.atac.filtered) <- "ACTIVITY"
# 
# markers <- FindAllMarkers(expr.atac.filtered, only.pos = T, min.pct = 0.1, logfc.threshold = 0.05)
# table(markers$cluster)
# markers_ftd <- markers[markers$p_val_adj <= 0.05,]
# table(markers_ftd$cluster)
# 
# top1 <- markers_ftd %>% group_by(cluster) %>% top_n(1, avg_logFC)
# top5 <- markers_ftd %>% group_by(cluster) %>% top_n(5, avg_logFC)
# 
# for(i in top1$gene) {
#     print(FeaturePlot(object = expr.atac.filtered, features = i, reduction = "tsne",
#                       pt.size = 0.5,label = T))
# }
# 
# write.table(markers_ftd,"Seurat_integration/cluster_marker_gene.txt",
#             row.names = T,col.names = T,sep='\t',quote=F)
# 
# DoHeatmap(expr.atac.filtered, features = top5$gene,disp.min = -3,disp.max = 3,size = 2.5)
# 
# dev.off()
