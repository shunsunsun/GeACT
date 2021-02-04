suppressMessages({
library(ArchR)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(Seurat)
library(plyr)
library(networkD3)
})

my_clip <- function(x, lb, ub){
  pmax(lb, pmin(x, ub))
}

change_assay <- function(sobj, assayin = "RNA", assayout = "ACTIVITY"){
  sobj[[assayout]] <- sobj[[assayin]]
  DefaultAssay(sobj) <- assayout
  sobj[[assayin]] <- NULL
  
  return(sobj)
}

coembedding <- function(sobj1, sobj2, features = NULL, disc = "tech", renorm = F,
                        assay1 = "RNA", assay2 = "RNA", subset = F, name1 = NULL, name2 = NULL){
  name1 <- if(is.null(name1)) assay1 else name1
  name2 <- if(is.null(name2)) assay2 else name2
  if(assay1 != "RNA") sobj1 <- change_assay(sobj1, assay1, "RNA")
  if(assay2 != "RNA") sobj2 <- change_assay(sobj2, assay2, "RNA")
  if (subset) {
    genes_both <- intersect(rownames(sobj1), rownames(sobj2))
    sobj1 <- sobj1[genes_both, ]
    sobj2 <- sobj2[genes_both, ]
  }
  merge_obj <- merge(sobj1, sobj2)
  cat(sprintf("Merge obj of (%d, %d) and (%d, %d) into (%d, %d)\n", 
              dim(sobj1)[1], dim(sobj1)[2],
              dim(sobj2)[1], dim(sobj2)[2],
              dim(merge_obj)[1], dim(merge_obj)[2]))
  
  if (renorm) merge_obj <- NormalizeData(merge_obj, normalization.method = "LogNormalize")
  features <- if (is.null(features)) VariableFeatures(merge_obj) else features
  merge_obj <- ScaleData(merge_obj, features = features)
  merge_obj <- RunPCA(merge_obj, features = features, npcs = 50, verbose = F)
  merge_obj <- RunTSNE(merge_obj, dims = 1:30)
  DimPlot(merge_obj, group.by = disc) + ggtitle(sprintf("Coembedding of %s and %s", name1, name2))
}

summarizedExperiment2Seurat <- function(se, assay = "RNA", project = "Project", do_norm = F){
  gene_meta <- as.data.frame(elementMetadata(se)) %>% column_to_rownames("name")
  gene_meta <- gene_meta[colnames(gene_meta) != "idx"]
  rownames(se) <- rownames(gene_meta)
  
  cell_meta <- as.data.frame(colData(se))
  
  se_assay <- assay(se)
  
  sobj <- CreateSeuratObject(counts = se_assay, project = project, meta.data = cell_meta, assay = assay)
  sobj[[assay]] <- AddMetaData(sobj[[assay]], gene_meta)
  
  if(do_norm) sobj <- NormalizeData(sobj, normalization.method = "LogNormalize")
  return(sobj)
}

integrate <- function(activity_atac, expr_rna, peak_matrix, activity_slot = "ACTIVITY", 
                      weight.reduction = NULL){
  # atac object setup
  activity_atac[["ATAC"]] <- CreateAssayObject(counts = assay(peak_matrix))
  activity_atac$tech <- "atac"
  
  DefaultAssay(activity_atac) <- activity_slot
  activity_atac <- NormalizeData(activity_atac)
  activity_atac <- ScaleData(activity_atac, features = rownames(activity_atac))
  
  DefaultAssay(activity_atac) <- "ATAC"
  VariableFeatures(activity_atac) <- names(which(Matrix::rowSums(activity_atac) > ncol(activity_atac)))
  activity_atac <- RunLSI(activity_atac, n = 50, scale.max = NULL)
  activity_atac <- RunUMAP(activity_atac, reduction = "lsi", dims = 1:50)
  # expr_atac <- NormalizeData(expr_atac)
  # expr_atac <- ScaleData(expr_atac, features = rownames(expr_atac[['ATAC']]))
  
  # rne object setup
  expr_rna <- NormalizeData(expr_rna)
  expr_rna[['ident']] <- Idents(expr_rna)
  expr_rna$tech <- "rna"
  genes.use <- intersect(rownames(expr_rna), rownames(activity_atac[[activity_slot]]))
  expr_rna_sub <- expr_rna[genes.use, ]
  expr_rna_sub <- FindVariableFeatures(expr_rna_sub, nfeatures = 2000)
  
  # find anchors and transfer data
  transfer.anchors <- FindTransferAnchors(
    reference = expr_rna_sub, query = activity_atac,
    features = VariableFeatures(object = expr_rna_sub),  #TOTEST 测试是否需要intersect
    normalization.method = "LogNormalize",
    reference.assay = "RNA", query.assay = activity_slot, reduction = "cca")
  
  if (is.null(weight.reduction)) weight.reduction <- activity_atac[["lsi"]]
  
  celltype.predictions <- TransferData(
    anchorset = transfer.anchors, refdata = expr_rna_sub$ident,
    weight.reduction = weight.reduction, dims = 1:length(weight.reduction))
  
  activity_atac <- AddMetaData(activity_atac, metadata = celltype.predictions$predicted.id, "predicted.id")
  activity_atac <- AddMetaData(activity_atac, metadata = celltype.predictions$prediction.score.max, "prediction.score.max")
  
  # to make the colors match
  activity_atac$predicted.id <- factor(activity_atac$predicted.id, levels = levels(expr_rna_sub))
  Idents(activity_atac) <- celltype.predictions$predicted.id
  # write.table(celltype.predictions,"../Seurat_integration/celltype_predictions.txt",
  #             sep='\t',quote = F,row.names = T,col.names = T)
  
  # Co-embedding
  refdata <- GetAssayData(expr_rna_sub, assay = "RNA", slot = "data")
  
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                             weight.reduction = weight.reduction, dims = 1:length(weight.reduction))
  cat(paste0("Tranfer dim: ", dim(imputation)[1], " ", dim(imputation)[2], "\n"))
  
  activity_atac[["RNA"]] <- imputation
  # saveRDS(expr.atac, file = "../Seurat_integration/Seurat_expr_atac_integration.rds")
  coembed <- merge(x = expr_rna_sub, y = activity_atac)
  
  # Finally, we run PCA and TSNE on this combined object, to visualize the co-embedding of both
  # datasets
  coembed <- ScaleData(coembed, features = genes.use, do.scale = F, verbose = F)
  coembed <- RunPCA(coembed, features = genes.use, verbose = F)
  coembed <- RunTSNE(coembed, dims = 1:30)
  coembed$celltype <- ifelse(!is.na(coembed$ident), coembed$ident, coembed$predicted.id)
  # saveRDS(coembed, file = "../Seurat_integration/Seurat_expr_coembed.rds")
  
  p1 <- DimPlot(coembed, group.by = "tech")
  p2 <- DimPlot(coembed, group.by = "ident", label = T, repel = T)
  print(p1 + p2 + ggtitle("Coembedding of integrated RNA and ATAC"))
  
  list(imputation = imputation, celltype.predictions = celltype.predictions)
}

peakMatrix2ciceroCDS <- function(peakMatrix){
  indata <- assay(peakMatrix)
  # binarize the matrix
  indata@x[indata@x > 0] <- 1
  
  # format cell info
  cellinfo <- as.data.frame(colData(peakMatrix))
  
  # format peak info
  peaks <- rowRanges(peakMatrix)
  peakinfo <- data.frame(row.names = paste(as.character(seqnames(peaks)), start(peaks), end(peaks), sep = "_"))
  
  rownames(peakMatrix) <- rownames(peakinfo)
  
  rownames(indata) <- rownames(peakinfo)
  colnames(indata) <- rownames(cellinfo)
  
  # make CDS
  input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                  cell_metadata = cellinfo,
                                                  gene_metadata = peakinfo))
  
  #Ensure there are no peaks included with zero reads
  dim(input_cds)
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]
  dim(input_cds)
  
  #### dimensionality reduction by UMAP
  set.seed(2020)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
  plot_cells(input_cds)
  
  #### access the UMAP coordinates
  umap_coords <- reducedDims(input_cds)$UMAP
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
  
  return(list(input_cds=input_cds, 
              cicero_cds=cicero_cds))
}

cmp_transfer <- function(transfer1, transfer2, lev, name1 = "geneScore", name2 = "cicero"){
  cat(sprintf("%d same prediction from %s and %s out of %d\n", 
              sum(transfer1 == transfer2), 
              name1, name2, 
              length(transfer1)))
  
  lev1 <- paste(name1, lev, sep = "_")
  lev2 <- paste(name2, lev, sep = "_")
  
  transfer1 <- factor(paste(name1, transfer1, sep = "_"), levels = lev1)
  transfer2 <- factor(paste(name2, transfer2, sep = "_"), levels = lev2)
  
  data <- data.frame(transfer1 = transfer1, transfer2 = transfer2)
  links <- ddply(data, c("transfer1", "transfer2"), nrow)
  colnames(links) <- c("source", "target", "values")
  
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(levels(links$source), 
           levels(links$target)), 
    group=rep(lev, 2)
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     NodeID = "name", Value = "values",
                     sinksRight=FALSE, NodeGroup = "group")
}

# cell group
ident2clgrp <- function(ident.input) {
  ident.input[grepl("^Epi", ident.input)] <- "Epithelial"
  ident.input[grepl("^Endo", ident.input)] <- "Endothelial"
  ident.input[grepl("^SM-", ident.input)] <- "Smooth muscle"
  ident.input[grepl("^SM$", ident.input)] <- "Smooth muscle"
  ident.input[grepl("^SKM$", ident.input)] <- "Skeletal muscle"
  ident.input[grepl("^Fibro", ident.input)] <- "Fibroblast"
  # immune
  ident.input[grepl("^B-", ident.input)] <- "B" 
  ident.input[grepl("^Pro-B", ident.input)] <- "B" 
  ident.input[grepl("^Pre-B", ident.input)] <- "B" 
  ident.input[grepl("^DC/Macro", ident.input)] <- "DC/Macrophage"
  ident.input[grepl("^Mast-", ident.input)] <- "Mast"
  ident.input[grepl("^Neutrophil-", ident.input)] <- "Neutrophil"
  ident.input[grepl("^NKT-", ident.input)] <- "NKT"
  ident.input[grepl("^T-", ident.input)] <- "T" 
  ident.input[grepl("^Pre-T", ident.input)] <- "T" 
  #
  ident.input[grepl("^Erythrocyte-", ident.input)] <- "Erythrocyte"
  ident.input[grepl("^CACNA1A-", ident.input)] <- "CACNA1A"
  ident.input[ident.input %in% c("PT", "LoH", "LoH-Prog", "DT", "PC-CLU", "PC-BCAT1", "Podocyte-GPC3", "Podocyte-PLA2R1")] <- "Epithelial"
  ident.input[grepl("^Sertoli-", ident.input)] <- "Sertoli"
  ident.input[grepl("^Granulosa-", ident.input)] <- "Granulosa"
  # FGC
  ident.input[grepl("^SSC$", ident.input)] <- "FGC"
  return(ident.input)
}

plotMarker <- function(proj, geneScore, genes, do_plot = T){
  figs <- list()
  for (gene in genes){
    if (!(gene %in% rownames(geneScore))) next
    f1 <- plotEmbedding(proj, embedding = "peakUMAP", colorBy = "GeneScoreMatrix", name = gene, plotAs = "points", size = 1.5, continuousSet = "whiteBlue", imputeWeights = NULL)
    f2 <- VlnPlot(geneScore, features = gene, group.by = "group")
    if (do_plot){
      print(f1)
      print(f2)
    }
    figs[[gene]] <- list(f1, f2)
  }
  figs
}

simpleMarkersPlot <- function(proj, geneScore){
  pdf("markers.pdf", width = 10, height = 10)
  print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
  
  # Epithelial
  plotMarker(proj, geneScore, "EPCAM")
  
  # Endothelial
  plotMarker(proj, geneScore, "PECAM1")
  
  # Smooth muscle
  plotMarker(proj, geneScore, "ACTG2")
  
  # Fibroblast
  plotMarker(proj, geneScore, c("COL1A1", "COL2A1"))
  
  # Glial
  plotMarker(proj, geneScore, "PLP1")
  
  # Immune
  plotMarker(proj, geneScore, "PTPRC")
  
  # T
  plotMarker(proj, geneScore, "CD3D")
  # B
  plotMarker(proj, geneScore, c("CD19", "CD79A"))
  
  # Erythrocyte
  plotMarker(proj, geneScore, "HBG1")
  
  # Proliferative
  plotMarker(proj, geneScore, "TOP2A")
  
  dev.off()
}


fullMarkersPlot <- function(proj, geneScore){
  dir.create("markers_plot", showWarnings = F)
  
  pdf("Epithelial.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, c("EPCAM", "VIM", "PTPRC", "HBG1", "PECAM1", "COL1A1"))
  dev.off()
  # epi, emt基质, immune, ery, endo,  
  
  pdf("Endothelial.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, "EPCAM")
  dev.off()
  
  pdf("Smooth_muscle.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, c("ACTG2", "CNN1", "ACTA2", "TAGLN", "NOTCH3", "LGR5", "DES", "LGR6", "MYH11"))
  dev.off()
  
  pdf("Glial.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, "EPCAM")
  dev.off()
  
  pdf("Fibroblast.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, c("COL1A1", "COL2A1", "PDGFRA", "ELN", "ACTA2", "PLIN2", "APOE"))
  dev.off()
  
  pdf("B.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, c("CD79A", "CD24", "MS4A1", "CD19"))
  dev.off()
  
  pdf("T.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, c("CD3D", "CCR7", "LEF1"))
  dev.off()
  
  pdf("DCMC.pdf", width = 10, height = 10)
  plotMarker(proj, geneScore, c("HLA-DRB1", "CLEC9A", "LAMP3", "CD1C", "PLD4", "CD14", "S100A8", "FCGR3A"))
  dev.off()
}
