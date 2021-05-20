suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(plyr)
  library(networkD3)
  library(GenomicRanges)
  library(assertthat)
  
  library(miniUI)
  library(shiny)
  
  library(RColorBrewer)
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

summarizedExperiment2Seurat <- function(se, assay = "RNA", project = "Project", do_norm = F, rename_assay = NULL){
  gene_meta <- as.data.frame(elementMetadata(se)) %>% column_to_rownames("name")
  gene_meta <- gene_meta[colnames(gene_meta) != "idx"]
  rownames(se) <- rownames(gene_meta)
  
  cell_meta <- as.data.frame(colData(se))
  
  se_assay <- assays(se)[[assay]]
  
  rename_assay <- if( is.null(rename_assay) ) assay else rename_assay
  
  sobj <- CreateSeuratObject(counts = se_assay, project = project, meta.data = if(ncol(cell_meta) > 0) cell_meta else NULL, 
                             assay = rename_assay)
  if(ncol(gene_meta) > 0) sobj[[rename_assay]] <- AddMetaData(sobj[[rename_assay]], gene_meta)
  
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
  
  tsne_data <- FetchData(coembed, vars = c("tSNE_1", "tSNE_2", "tech", "celltype"))
  
  p1 <- DimPlot(coembed, group.by = "tech")
  p2 <- DimPlot(coembed, group.by = "ident", label = T, repel = T)
  print(p1 + p2 + ggtitle("Coembedding of integrated RNA and ATAC"))
  
  list(imputation = imputation, celltype.predictions = celltype.predictions, tsne.data = tsne_data)
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
  ident.input[grepl("^Macro", ident.input)] <- "Macrophage"
  ident.input[grepl("^Neutrophil-", ident.input)] <- "Neutrophil"
  ident.input[grepl("^NKT-", ident.input)] <- "NKT"
  ident.input[grepl("^T-", ident.input)] <- "T" 
  ident.input[grepl("^Pre-T", ident.input)] <- "T" 
  ident.input[grepl("^Immune-", ident.input)] <- "Immune"
  #
  ident.input[grepl("^Erythrocyte-", ident.input)] <- "Erythrocyte"
  ident.input[grepl("^CACNA1A-", ident.input)] <- "CACNA1A"
  ident.input[ident.input %in% c("PT", "LoH", "LoH-Prog", "DT", "PC-CLU", "PC-BCAT1", "Podocyte-GPC3", "Podocyte-PLA2R1")] <- "Epithelial"
  ident.input[ident.input %in% c("LSEC")] <- "Endothelial"
  ident.input[grepl("^Sertoli-", ident.input)] <- "Sertoli"
  ident.input[grepl("^Granulosa-", ident.input)] <- "Granulosa"
  # FGC
  ident.input[grepl("^SSC$", ident.input)] <- "FGC"
  return(ident.input)
}

# standardize genes to deal with seurat substitution of "_" with "_" problem
standardize_genename <- function(genes, genes_nomenclature) {
  genes_nomenclature[match(gsub("_", "-", genes), gsub("_", "-", genes_nomenclature))]
}

plotMarker <- function(proj, genes, geneScore = NULL, do_plot = T, group_by = "group"){
  figs <- list()
  
  if( is.null(geneScore) ) geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = F) %>% 
      summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")
  
  # check genes to deal with seurat substitution of "_" with "_" problem
  genes_anno <- getGeneAnnotation(proj)$genes$symbol
  genes <- standardize_genename(genes, genes_anno)
  
  for (gene in genes){
    if (!any(c(gene, gsub(gene, "_", "-")) %in% rownames(geneScore))) next
    
    # trycatch block to deal with seurat substitution of "_" with "_" problem
    f1 <- plotEmbedding(proj, embedding = "peakUMAP", colorBy = "GeneScoreMatrix", name = gene, plotAs = "points", size = 1.5, continuousSet = "whiteBlue", imputeWeights = NULL)
    
    f2 <- tryCatch({
      VlnPlot(geneScore, features = gene, group.by = group_by)
    }, 
    error = function(cond) {
      VlnPlot(geneScore, features = gsub("_", "-", gene), group.by = group_by)
    }
    )
    
    # f3 <- plotBrowserTrack(
    #   ArchRProj = proj, 
    #   groupBy = group_by, 
    #   geneSymbol = gene, 
    #   upstream = 100000,
    #   downstream = 100000,
    #   baseSize = 12,
    #   facetbaseSize = 12,
    #   loops = NULL # getCoAccessibility(proj_epi) # getPeak2GeneLinks(proj_epi)
    # )
    if (do_plot){
      print(f1)
      print(f2)
      # grid::grid.newpage()
      # grid::grid.draw(f3[[1]])
    }
    figs[[gene]] <- list(f1, f2) #, f3)
  }
  figs
}

# plot general markers
simpleMarkersPlot1 <- function(proj, geneScore = NULL){
  print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
  
  # for performance issue，take geneScore as a parameter
  if( is.null(geneScore) ) geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = F) %>% 
      summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")
  
  # Epithelial
  grid::grid.newpage()
  grid::grid.text("Epithelial", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("EPCAM", "KRT5"), geneScore)
  
  # Endothelial
  grid::grid.newpage()
  grid::grid.text("Endothelial", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("PECAM1", "VWF", "CDH13"), geneScore)
  
  # Smooth muscle
  grid::grid.newpage()
  grid::grid.text("Smooth muscle", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("ACTG2", "ACTA2", "CNN1", "TAGLN", "MYH11"), geneScore)
  
  # Fibroblast
  grid::grid.newpage()
  grid::grid.text("Fibroblast", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("COL1A1", "COL2A1", "PDGFRA"), geneScore)
  
  # Pericyte
  grid::grid.newpage()
  grid::grid.text("Pericyte", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("CSPG4", "TRPC6", "PDGFRB"), geneScore)
  
  # Glial
  grid::grid.newpage()
  grid::grid.text("Glial", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("PLP1", "S100B"), geneScore)
  
  # Immune
  grid::grid.newpage()
  grid::grid.text("Immune", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, "PTPRC", geneScore)
  
  # T
  grid::grid.newpage()
  grid::grid.text("T", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, "CD3D", geneScore)
  
  # NKT
  grid::grid.newpage()
  grid::grid.text("NKT", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("CD3E", "CD8A", "FCER1G", "TYROBP"), geneScore)
  
  # B
  grid::grid.newpage()
  grid::grid.text("B", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("CD19", "CD79A", "CD24", "MS4A1"), geneScore)
  
  # DC
  # grid::grid.newpage()
  # grid::grid.text("DC", gp = gpar(fontsize = 32, fontface = "bold"))
  # plotMarker(proj, "MHCII", geneScore)
  
  # Macrophage
  grid::grid.newpage()
  grid::grid.text("Macrophage", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("MARCO", "MSR1", "MRC1"), geneScore)
  
  # NK
  grid::grid.newpage()
  grid::grid.text("NK", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("KLRD1", "NKG7", "TYROBP"), geneScore)
  
  # Mast
  grid::grid.newpage()
  grid::grid.text("Mast", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("MS4A2", "CPA3", "TPSAB1"), geneScore)
  
  # Erythrocyte
  grid::grid.newpage()
  grid::grid.text("Erythrocyte", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("HBB", "HBG1", "HBG2"), geneScore)
  
  # Proliferative
  grid::grid.newpage()
  grid::grid.text("Proliferative", gp = grid::gpar(fontsize = 32, fontface = "bold"))
  plotMarker(proj, c("TOP2A", "MKI67"), geneScore)
  
  return(NULL)
}

# plot general markers to compare with RNA dotplot
simpleMarkersPlot2 <- function(proj, geneScore = NULL, tissue = "01_stomach", root = "/data/Lab/otherwork/GeACT/ATAC"){
  print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
  
  # for performance issue，take geneScore as a parameter
  if( is.null(geneScore) ) {
    geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = F) %>% 
      summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")
  }
  
  tissue_marker_table <- read.table( paste(root, "meta/marker_genes/marker_genes_ploted.txt", sep = "/"), sep = "\t", header = F, stringsAsFactors = F)
  
  markers_to_plot <- tissue_marker_table %>% filter(V1 == tissue) %>% .[[2]]
  
  cat(markers_to_plot, "\n")
  plotMarker(proj, markers_to_plot, geneScore)
  
  Idents(geneScore) <- geneScore[["predictedIdent"]]
  print(DotPlot(geneScore, features = markers_to_plot, idents = ) + RotatedAxis())
  Idents(geneScore) <- geneScore[["group"]]
  print(DotPlot(geneScore, features = markers_to_plot, idents = ) + RotatedAxis())
  
  return(NULL)
}

# plot signature/marker genes from DE analysis of RNA
signatureGenesPlot <- function(proj, geneScore = NULL, tissue = "01_stomach", root = "/data/Lab/otherwork/GeACT/ATAC") {
  print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "group", plotAs = "points", size = 1.5))
  
  # for performance issue，take geneScore as a parameter
  if( is.null(geneScore) ) geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = F) %>% 
      summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")
  
  tissue_signature_table <- read.table( paste(root, "meta/marker_genes/cellType_signatureGenes.txt", sep = "/"), sep = "\t", header = T, stringsAsFactors = F)
  
  tissue_raw <- gsub("^.*?_", "", tissue) %>% gsub("_", " ", .)
  signatures_df <- tissue_signature_table %>% filter(tissue == tissue_raw) %>% group_by(ident) %>% top_n(n = 2, wt = myAUC) %>% select(ident, gene)
  
  print(signatures_df)
  for(.ident in unique( signatures_df$ident )){
    signature_genes <- signatures_df %>% subset(ident == .ident) %>% .$gene
    
    grid::grid.newpage()
    grid::grid.text(.ident, gp = grid::gpar(fontsize = 32, fontface = "bold"))
    plotMarker(proj, signature_genes, geneScore)
  }
  
  return(NULL)
}

tunedSignatureGenesPlot <- function(proj, topn_df, geneScore = NULL) {
  print(plotEmbedding(proj, embedding = "peakUMAP", colorBy = "cellColData", name = "tuned_ident", plotAs = "points", size = 1.5))
  
  # for performance issue，take geneScore as a parameter
  if( is.null(geneScore) ) geneScore <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = F) %>% 
      summarizedExperiment2Seurat(assay = "GeneScoreMatrix", rename_assay = "ACTIVITY")
  
  for(.cluster in levels(topn_df$cluster)){
    cluster_signatures <- topn_df %>% filter(cluster == .cluster) %>% pull(gene)
    
    grid::grid.newpage()
    grid::grid.text(paste0("cluster: ", .cluster), gp = grid::gpar(fontsize = 32, fontface = "bold"))
    plotMarker(proj, genes = cluster_signatures, geneScore = geneScore, group_by = "tuned_ident")
  }
}

extend_range <- function(lower_bound, upper_bound, extend_ratio = 0.1) {
  range_len <- upper_bound - lower_bound
  c(lower_bound - extend_ratio * range_len, upper_bound + extend_ratio * range_len)
}

# helper function, data_df contains x, y coordination, and features to plot
# feature should be a continuous value column, 
# group_bys should be discrete/factor-like value columns
interactive_choose_df <- function(data_df, features = NULL, point_size = 1.5, group_bys = NULL) {
  assert_that(all(c("x", "y") %in% colnames(data_df)))
  
  if (!is.null(group_bys)) {
    assert_that(all(group_bys %in% colnames(data_df)))
  }
  
  if (!is.null(features)) {
    assert_that(all(features %in% colnames(data_df)))
  }
  
  need_choose_color_by <- (!is.null(features)) | (!is.null(group_bys))
  if(need_choose_color_by) {
    color_by_options <- c("none", if( !is.null(features) ) "features" else NULL,
                          if( !is.null(group_bys) ) "groups" else NULL)
  }
  
  # helper function for color scales
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  
  # generage color scales
  color_scales <- lapply(group_bys, function(.group_by) {
    num_of_colors <- data_df[[.group_by]] %>% unique %>% length
    color_scale <- getPalette(num_of_colors)
    names(color_scale) <- data_df[[.group_by]] %>% unique
    return(color_scale)
  })
  names(color_scales) <- group_bys
  
  ui <- miniPage(
    gadgetTitleBar("Select Cells"), 
    fillRow(flex = c(5, 1), 
            miniContentPanel(
              padding = 0,
              plotOutput("plot1", height = "100%", brush = "brush", click = "plot_click")
            ),
            miniContentPanel(
              if (need_choose_color_by) selectInput("color_by", "Color by features, groups, or none", choices = color_by_options) else NULL, 
              
              if (!is.null(features)) uiOutput("chooseFeatures") else NULL, 
              
              if (!is.null(group_bys)) selectInput("group_by", "Group cells by", choices = group_bys) else NULL, 
              if (!is.null(group_bys)) uiOutput("visableGroups") else NULL, 
              
              htmlOutput('click_info'),
              htmlOutput('brush_info')
            )
    ), 
    miniButtonBlock(
      actionButton("add", "", icon = icon("thumbs-up")),
      actionButton("sub", "", icon = icon("thumbs-down")),
      actionButton("none", "", icon = icon("refresh")),
      # actionButton("all", "", icon = icon("refresh"))
    )
  )
  
  server <- function(input, output) {
    # TODO, when default values are not selected, how to do plot problem, including input$visable_group, input$feature, etc. any better realizations than req implementation?
    
    # dynamic UIs
    output$visableGroups <- renderUI({
      groups_to_choose <- c("all", data_df[[input$group_by]] %>% unique)
      selectInput("visable_group", "Groups to visualize", choices = groups_to_choose, multiple = T, selected = "all")
    })
    
    output$chooseFeatures <- renderUI(({
      # need_choose_color_by is true gurantees that input$color_by exists
      if(need_choose_color_by && input$color_by == "features") selectInput("feature", "Feature to colorby", 
                                                                           choices = features) else NULL
    }))
    
    
    # For storing selected points
    vals <- reactiveValues(keep = rep(F, nrow(data_df)), xy = c(0, 0))
    
    observeEvent(input$done, stopApp(vals$keep))
    
    # compute x, y range for plots, in case it changes according to subsets
    x_range <- extend_range(data_df$x %>% min, data_df$x %>% max)
    y_range <- extend_range(data_df$y %>% min, data_df$y %>% max)
    
    # remove all legends for consistent visualization
    
    output$plot1 <- renderPlot({
      # Plot the kept and excluded points as two separate data sets
      keep <- data_df[vals$keep, , drop = F]
      exclude <- data_df[!vals$keep, , drop = F]
      
      if(!is.null(group_bys)) {
        req(input$visable_group)
        
        if ("all" %in% input$visable_group){
        } else {
          keep    <- keep[keep[[input$group_by]] %in% input$visable_group, , drop = F]
          exclude <- exclude[exclude[[input$group_by]] %in% input$visable_group, , drop = F]
        }
      }
      
      fig <- ggplot(exclude, aes(x, y)) + 
        xlim(x_range[1], x_range[2]) + 
        ylim(y_range[1], y_range[2]) + 
        theme(
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          
          panel.background = element_rect(fill = 'white', colour = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          axis.line = element_blank(),
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          
          legend.position = "none"
        ) + 
        geom_point(data = keep, color = "green", size = point_size)
      
      if(need_choose_color_by) req(input$color_by)
      if(need_choose_color_by & input$color_by == "features") {
        req(input$feature)
        fig <- fig + geom_point(size = point_size, mapping = aes_string(colour = input$feature)) + 
          scale_color_gradient(low = "#e7e8ea", high = "blue", limits = c(0, max(data_df[[input$feature]]) ) ) + 
          ggtitle(input$feature)
      } else if (need_choose_color_by & input$color_by == "groups") {
        req(input$group_by)
        fig <- fig + geom_point(size = point_size, mapping = aes_string(colour = input$group_by)) + 
          scale_color_manual(values = color_scales[[input$group_by]]) + 
          ggtitle(input$group_by)
      } else if (!need_choose_color_by || input$color_by == "none") {
        fig <- fig + geom_point(size = point_size) + 
          ggtitle("default")
      } else {
        stopApp(stop("Oops, something goes wrong!!!"))
      }
      
      fig
    })
    
    # Update selected points
    selected <- reactive({
      if(!is.null(input$visable_group) & !("all" %in% input$visable_group) ) brushedPoints(data_df, input$brush, allRows = T)$selected_ & (data_df[[input$group_by]] %in% input$visable_group ) else 
        brushedPoints(data_df, input$brush, allRows = T)$selected_
    })
    
    observeEvent(input$add, vals$keep <- vals$keep | selected())
    observeEvent(input$sub, vals$keep <- vals$keep & !selected())
    observeEvent(input$none, vals$keep <- rep(F, nrow(data_df)))
    # observeEvent(input$all, vals$keep <- rep(T, nrow(data_df)))
    
    # coordinate info
    observeEvent(input$plot_click, vals$xy <- c(input$plot_click$x, input$plot_click$y))
    
    output$click_info <- renderUI({
      paste0(sprintf("x: %.2f<br/>", vals$xy[1]),
             sprintf("y: %.2f", vals$xy[2]) ) %>% HTML
    })
    
    output$brush_info <- renderUI({
      req(input$brush)
      paste0(sprintf("x: [%.2f, %.2f]<br/>", input$brush$xmin, input$brush$xmax), 
             sprintf("y: [%.2f, %.2f]",   input$brush$ymin, input$brush$ymax)) %>% HTML
    })
  }
  
  runGadget(ui, server, viewer = dialogViewer("choose_cell", width = 1000, height = 1000))
}

# bridge function between ArchR proj and interactive_choose_cell
interactive_choose_cell <- function(proj, features = NULL, embedding = "peakUMAP", log2Norm = F, group_bys = NULL) {
  data_df <- getEmbedding(proj, embedding = embedding)
  names(data_df) <- c("x", "y")
  if(!is.null(feature)) {
    feature_df <- ArchR:::.getMatrixValues(proj, name = features, matrixName = "GeneScoreMatrix") %>% t %>% as.data.frame
    if(log2Norm) feature_df <- log2(feature_df + 1)
    data_df <- merge(data_df, feature_df, by = 0) %>% column_to_rownames("Row.names")
  }
  
  if(!is.null(group_bys)) {
    meta_df <- getCellColData(proj, select = group_bys)
    data_df <- merge(data_df, meta_df %>% as.data.frame, by = 0) %>% column_to_rownames("Row.names")
  }
  
  rownames(data_df[interactive_choose_df(data_df, features, group_bys = group_bys), , drop = F])
}

# choose cell based on embedding, static fine tuning tool for cell type annotation
static_choose_cell <- function(proj, embedding = "peakUMAP", xlimit = c(0, 10), ylimit = c(0, 10), group_by = NULL, used_groups = NULL) {
  assert_that(length(xlimit) == 2 && is.numeric(xlimit) && xlimit[1] <= xlimit[2])
  assert_that(length(ylimit) == 2 && is.numeric(ylimit) && ylimit[1] <= ylimit[2])
  
  data_df <- getEmbedding(proj, embedding = embedding)
  names(data_df) <- c("x", "y")
  
  if(!is.null(group_by)) {
    meta_df <- getCellColData(proj, select = group_by)
    data_df <- merge(data_df, meta_df %>% as.data.frame, by = 0) %>% column_to_rownames("Row.names")
    
    data_df <- data_df[data_df[[group_by]] %in% used_groups, , drop = T]
  }
  
  data_df[data_df$x >= xlimit[1] & data_df$x <= xlimit[2] &
           data_df$y >= ylimit[1] & data_df$y <= ylimit[2], , drop = F] %>% rownames
}

chromVARFeaturePlot <- function(object, feature, data = NULL, reduction = NULL, dims = c(1, 2), cells = NULL, size = 3) {
  if (is.null(data)){
    reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", feature), cells = cells, slot = "data")
  } else {
    dims <- colnames(data)
    reduction <- "UMAP"
    feature_df <- FetchData(object = object, vars = feature, cells = cells, slot = "data")
    data <- merge(data, feature_df, by.x = 0, by.y = 0)
  }
  
  ggplot(data = data, mapping = aes_string(x = paste0("`", dims[1], "`"), y = paste0("`", dims[2], "`"), colour = paste0("`", feature, "`"))) + geom_point(size = 1) +
    scale_color_gradientn(colours = c("blue","#f1f1f1","red"), 
                          values = scales::rescale(c(min(data[[feature]]), 0, max(data[[feature]]))), 
                          guide = "colorbar",
                          limits = c(min(data[[feature]]), max(data[[feature]]))) +
    ggtitle(feature) +# ggtitle(gsub("-", "_", feature)) +
    xlab(paste(toupper(reduction), "1", sep = "_")) + ylab(paste(toupper(reduction), "2", sep = "_")) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      
      panel.background = element_rect(fill = 'white', colour = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      
      legend.background = element_rect(fill = "white", size = size, colour = "white"),
      legend.title = element_blank(),
      #legend.justification = c(0, 1),
      axis.line = element_blank(), # element_line(colour = "black"),
      axis.title = element_blank(), 
      axis.ticks = element_blank(), # element_line(colour = "black", size = 0.2),
      axis.text = element_blank()
    )
}

my_genomicDensity <- function (region, window.size = 1e7, n.window = NULL, overlap = TRUE, 
                               count.by = "value", chr.len = NULL, return.log = FALSE, max.clip = NULL) {
  chrs <- unique(seqnames(region))
  if(length(chrs) > 1){
    return(do.call("rbind", lapply(chrs, function(chr){
      if (is.null(chr.len)) {
        max_rg = NULL
      } else {
        if (chr %in% names(chr.len)) {
          max_rg = chr.len[chr]
        } else {
          max_rg = NULL
        }
      }
      df = my_genomicDensity(region[seqnames(region) == chr, ], 
                             window.size = window.size, overlap = overlap, 
                             chr.len = max_rg, count.by = count.by, return.log = return.log, max.clip = max.clip)
      # cbind(chr = rep(chr, nrow(df)), df)
    })))
  }
  region <- GenomicRanges::sort(region)
  # region <- GenomicRanges::reduce(region)
  
  if (!is.null(chr.len)) {
    max_pos <- max(c(chr.len, max(end(region))))
  } else {
    max_pos <- max(end(region))
  }
  
  if (overlap) {
    if (missing(n.window)) {
      b = seq(1, max_pos, by = window.size/2)
      s = b[-length(b)]
      s = s[-length(s)]
      e = s + window.size - 1
    }
    else {
      b = seq(1, max_pos, length = 2 * n.window - 1)
      s = b[-length(b)]
      s = s[-length(s)]
      e = s + b[3] - b[1] - 1
    }
  }
  else {
    if (missing(n.window)) {
      b = seq(1, max_pos, by = window.size)
      s = b[-length(b)]
      e = s + window.size - 1
    }
    else {
      b = seq(1, max_pos, length = n.window)
      s = b[-length(b)]
      e = s + b[2] - b[1]
    }
  }
  
  s = as.integer(s)
  e = as.integer(e)
  y = rep(0, length(s))
  
  names(y) = paste(s, e, sep = ",")
  windows = data.frame(start = s, end = e)
  
  region <- as.data.frame(region)[, c("start", "end", count.by)]
  op <- my_overlap_region(windows, region, count_by = count.by)
  
  if(return.log) op <- log(op + 1)
  if (!is.null(max.clip)) op <- my_clip(op, 0, max.clip)
  res <- data.frame(chr = chrs, start = s, end = e, value = op)
  return(res)
  
}

my_overlap_region <- function (gr1, gr2, count_by = "value") {
  nr1 = nrow(gr1)
  nr2 = nrow(gr2)
  overlap = rep(0, length = nr1)
  if (nr1 == 0) {
    return(overlap)
  }
  k_gr2 = 1
  for (i in seq_len(nr1)) {
    for (j in seq(k_gr2, nr2)) {
      if (gr2[j, 2] < gr1[i, 1]) {
        k_gr2 = ifelse(k_gr2 < nr2, k_gr2 + 1, nr2)
        next
      }
      else if (gr2[j, 1] > gr1[i, 2]) {
        break
      }
      else {
        # cat(sprintf("j: %d\n", j))
        # print(gr2[j, , drop = F])
        overlap[i] = overlap[i] + gr2[j, count_by]
      }
    }
  }
  return(overlap)
}

df2gr <- function(df){
  gr <- makeGRangesFromDataFrame(setNames(df, c("chr", "start", "end")))
  mcols(gr) <- df[, -(1:3), drop = F]
  
  return(gr)
}

gr2df <- function(gr, keep_strand_width = F) {
  df <- as.data.frame(gr)
  if (!keep_strand_width) df[, c("width", "strand")] <- NULL
  colnames(df)[1] <- "chr"
  
  return(df)
}

plotCircosFromRangeSE <- function(peak_matrix, group_by = "group", groups_to_cmp = NULL, window.size = 1e6, max.clip = 30, 
                                  cytoband = NULL, chrs = 1:22, use_colors = NULL, genes_df = NULL){
  if (is.null(groups_to_cmp)) groups_to_cmp  <- colData(peak_matrix)[, group_by] %>% table %>% sort(decreasing = T) %>% head(7) %>% names
  
  pseu_bulk_gr_list <- lapply(groups_to_cmp, function(group){
    group_peak_matrix <- peak_matrix[, peak_matrix$group == group]
    group_pseu_bulk <- rowMeans(assay(group_peak_matrix))
    group_pseu_bulk_gr <- rowRanges(group_peak_matrix)
    group_pseu_bulk_gr$value <- group_pseu_bulk
    return(group_pseu_bulk_gr)
  })
  
  names(pseu_bulk_gr_list) <- groups_to_cmp
  
  if (!is.null(cytoband)){
    if (is.data.frame(cytoband)){
      cytoband_df <- cytoband
    } else {
      cytoband_df <- cytoband$df
    }
    pseu_bulk_gr_list <- lapply(pseu_bulk_gr_list, function(pseu_bulk_gr) subsetByOverlaps(pseu_bulk_gr, df2gr(cytoband_df)))
  }
  
  plot_dfs <- lapply(pseu_bulk_gr_list, my_genomicDensity, window.size = window.size, return.log = F, max.clip = max.clip) # 1e6 30 # 5e4 10
  if (!is.null(cytoband)) plot_dfs <- lapply(plot_dfs, function(plot_df) subsetByOverlaps(df2gr(plot_df), df2gr(cytoband_df)) %>% gr2df)
  
  
  if (is.null(use_colors)) use_colors <- brewer.pal(length(groups_to_cmp), name = "Set2")
  
  y_max <- sapply(plot_dfs, function(df) {max(df$value)}) %>% max
  
  circos.par(start.degree = 90, gap.degree = 3)
  
  if (is.null(cytoband)){
    circos.initializeWithIdeogram(species = "hg38", plotType = c("axis", "labels"), chromosome.index = paste0("chr", chrs))
  } else {
    circos.initializeWithIdeogram(cytoband = cytoband_df, plotType = c("axis", "labels"), tickLabelsStartFromZero = F, chromosome.index = paste0("chr", chrs))
  }
  
  
  if (!is.null(genes_df)){
    gene_colors <- brewer.pal(nrow(genes_df), name = "Pastel1")
    circos.genomicTrack(genes_df, track.height = 0.1, ylim = c(0, 1), bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          i = getI(...)
                          circos.genomicRect(region, value, border = NA, col = gene_colors[i], ...)
                          circos.genomicText(region, value, labels = value[["name"]], y = 0.5)
                        }
    )
  }
  
  for(i in seq_along(pseu_bulk_gr_list)){
    set_track_gap(mm_h(2))
    circos.genomicTrack(plot_dfs[[i]], track.height = 0.1, ylim = c(0, y_max), bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value, type = "l", border = NA, baseline = 0, area = T,
                                              col = use_colors[i])
                        })
  }
  
  set_track_gap(mm_h(5))
  
  circos.genomicIdeogram(track.height = 0.03)
  circos.clear()
  
  lgd_points <- Legend(at = names(use_colors), type = "points", 
                       legend_gp = gpar(col = use_colors), title_position = "topleft", 
                       title = "Cell Group")
  draw(lgd_points)
}

# modified version of addGeneIntegrationMatrix in ArchR, to support passing arguments
# to seurat cca transfer process
my_addGeneIntegrationMatrix <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = NULL,
  groupATAC = NULL,
  groupRNA = NULL,
  groupList = NULL,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  embeddingATAC = NULL,
  embeddingRNA = NULL,
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  plotUMAP = TRUE,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE),
  nGenes = 2000,
  useImputation = TRUE,
  reduction = "cca",
  addToArrow = TRUE,
  scaleTo = 10000,
  genesUse = NULL,
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  transferParams = list(),
  threads = getArchRThreads(),
  verbose = TRUE,
  force = FALSE,
  logFile = createLogFile("addGeneIntegrationMatrix"),
  k.weight = 50,
  ...
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment", "Seurat"))
  .validInput(input = groupATAC, name = "groupATAC", valid = c("character", "null"))
  .validInput(input = groupRNA, name = "groupRNA", valid = c("character"))
  .validInput(input = groupList, name = "groupList", valid = c("list", "null"))
  .validInput(input = sampleCellsATAC, name = "sampleCellsATAC", valid = c("integer", "null"))
  .validInput(input = sampleCellsRNA, name = "sampleCellsRNA", valid = c("integer", "null"))
  .validInput(input = embeddingATAC, name = "embeddingATAC", valid = c("data.frame", "null"))
  .validInput(input = embeddingRNA, name = "embeddingRNA", valid = c("data.frame", "null"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = plotUMAP, name = "plotUMAP", valid = c("boolean"))
  .validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))  
  .validInput(input = nGenes, name = "nGenes", valid = c("integer"))
  .validInput(input = useImputation, name = "useImputation", valid = c("boolean"))
  .validInput(input = reduction, name = "reduction", valid = c("character"))
  .validInput(input = addToArrow, name = "addToArrow", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = genesUse, name = "genesUse", valid = c("character", "null"))
  .validInput(input = nameCell, name = "nameCell", valid = c("character"))
  .validInput(input = nameGroup, name = "nameGroup", valid = c("character"))
  .validInput(input = nameScore, name = "nameScore", valid = c("character"))
  .validInput(input = transferParams, name = "transferParams", valid = c("list"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logDiffTime("Running Seurat's Integration Stuart* et al 2019", tstart, verbose = verbose, logFile = logFile)
  
  .requirePackage("Seurat", source = "cran")
  
  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "Input-Parameters", logFile=logFile)
  
  if(is.null(groupList)){ #If null use all cells (blocking will still occur)
    groupList <- SimpleList()
    groupList[[1]] <- SimpleList(
      ATAC = ArchRProj$cellNames,
      RNA = colnames(seRNA)
    )
  }
  
  #########################################################################################
  # 1. Check All ATAC is Accounted For!
  #########################################################################################
  .logDiffTime("Checking ATAC Input", tstart, verbose = verbose, logFile = logFile)
  
  if(!is.null(groupATAC)){
    dfATAC <- getCellColData(ArchRProj = ArchRProj, select = groupATAC, drop = FALSE)
  }
  nCell <- rep(0, length(ArchRProj$cellNames))
  names(nCell) <- ArchRProj$cellNames
  
  groupList <- lapply(seq_along(groupList), function(x){
    
    ATAC <- groupList[[x]]$ATAC
    
    if(!is.null(groupATAC)){
      
      if(any(ATAC %in% dfATAC[,1])){
        idx <- which(ATAC %in% dfATAC[,1])
        ATAC2 <- rownames(dfATAC)[which(dfATAC[,1] %in% ATAC[idx])]
        if(length(idx) == length(ATAC)){
          ATAC <- ATAC2
        }else{
          ATAC <- c(ATAC[-idx], ATAC2)
        }
      }
      
    }
    
    SimpleList(ATAC = ATAC, RNA = groupList[[x]]$RNA)
    
  }) %>% SimpleList
  
  for(i in seq_along(groupList)){
    nCell[groupList[[i]]$ATAC] <- nCell[groupList[[i]]$ATAC] + 1
  }
  
  if(!all(nCell == 1)){
    .logMessage(paste0("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!"), logFile = logFile)
    stop("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!")
  }
  
  #########################################################################################
  # 2. Check All RNA is a Cell Name 
  #########################################################################################
  .logDiffTime("Checking RNA Input", tstart, verbose = verbose, logFile = logFile)
  
  #Set up RNA
  if(inherits(seRNA, "SummarizedExperiment")){
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if(groupRNA %ni% colnames(colData(seRNA))){
      .logMessage("groupRNA not in colData of seRNA", logFile = logFile)
      stop("groupRNA not in colData of seRNA")
    }
    seuratRNA$Group <- paste0(colData(seRNA)[, groupRNA, drop = TRUE])
    rm(seRNA)
  }else{
    if(groupRNA %ni% colnames(seRNA@meta.data)){
      .logMessage("groupRNA not in meta.data of Seurat Object", logFile = logFile)
      stop("groupRNA not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- paste0(seRNA@meta.data[,groupRNA])
    rm(seRNA)
  }
  gc()
  
  if(!is.null(groupRNA)){
    dfRNA <- DataFrame(row.names = colnames(seuratRNA), Group = seuratRNA$Group)
  }
  
  groupList <- lapply(seq_along(groupList), function(x){
    
    RNA <- groupList[[x]]$RNA
    
    if(!is.null(groupRNA)){
      
      if(any(RNA %in% dfRNA[,1])){
        idx <- which(RNA %in% dfRNA[,1])
        RNA2 <- rownames(dfRNA)[which(dfRNA[,1] %in% RNA[idx])]
        if(length(idx) == length(RNA)){
          RNA <- RNA2
        }else{
          RNA <- c(RNA[-idx], RNA2)
        }
      }
      
    }
    
    SimpleList(ATAC = groupList[[x]]$ATAC, RNA = RNA)
    
  }) %>% SimpleList
  
  cellRNA <- unlist(lapply(groupList, function(x) x$RNA))
  if(!all(cellRNA %in% colnames(seuratRNA))){
    .logMessage("Found cells for RNA not in colnames(seRNA)! Please retry your input!", logFile = logFile)
    stop("Found cells for RNA not in colnames(seRNA)! Please retry your input!")
  }
  
  seuratRNA <- seuratRNA[, unique(cellRNA)]
  seuratRNA <- NormalizeData(object = seuratRNA, verbose = FALSE)
  
  #########################################################################################
  # 3. Create Integration Blocks
  #########################################################################################
  
  #Check Gene Names And Seurat RowNames
  geneDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  sumOverlap <- sum(unique(geneDF$name) %in% unique(rownames(seuratRNA)))
  if(sumOverlap < 5){
    stop("Error not enough overlaps (",sumOverlap,") between gene names from gene scores (ArchR) and rna matrix (seRNA)!")
  }
  .logDiffTime(paste0("Found ", sumOverlap, " overlapping gene names from gene scores and rna matrix!"), tstart, verbose = TRUE, logFile = logFile)
  
  .logDiffTime("Creating Integration Blocks", tstart, verbose = verbose, logFile = logFile)
  
  blockList <- SimpleList()
  
  for(i in seq_along(groupList)){
    
    gLi <- groupList[[i]]
    
    #######################################
    # ATAC
    #######################################
    
    if(length(gLi$ATAC) > sampleCellsATAC){
      
      if(!is.null(embeddingATAC)){
        probATAC <- .getDensity(embeddingATAC[gLi$ATAC,1], embeddingATAC[gLi$ATAC,2])$density
        probATAC <- probATAC / max(probATAC)
        cellsATAC <- gLi$ATAC[order(probATAC, decreasing = TRUE)]
      }else{
        cellsATAC <- sample(gLi$ATAC, length(gLi$ATAC))
      }
      
      cutoffs <- lapply(seq_len(1000), function(x) length(gLi$ATAC) / x) %>% unlist
      blockSize <- ceiling(min(cutoffs[order(abs(cutoffs - sampleCellsATAC))[1]] + 1, length(gLi$ATAC)))
      
      #Density Based Blocking
      nBlocks <- ceiling(length(gLi$ATAC) / blockSize)
      
      blocks <- lapply(seq_len(nBlocks), function(x){
        cellsATAC[seq(x, length(cellsATAC), nBlocks)]
      }) %>% SimpleList
      
    }else{
      
      blocks <- list(gLi$ATAC)
    }
    
    #######################################
    # RNA
    #######################################
    
    if(!is.null(embeddingRNA)){
      probRNA <- .getDensity(embeddingRNA[gLi$RNA,1], embeddingRNA[gLi$RNA,2])$density
      probRNA <- probRNA / max(probRNA)
    }else{
      probRNA <- rep(1, length(gLi$RNA))
    }
    
    blockListi <- lapply(seq_along(blocks), function(x){
      
      SimpleList(
        ATAC = blocks[[x]],
        RNA = sample(x = gLi$RNA, size = min(sampleCellsRNA, length(gLi$RNA)) , prob = probRNA)
      )
      
    }) %>% SimpleList
    
    blockList <- c(blockList, blockListi)
    
  }
  rm(groupList)
  
  #########################################################################################
  # 4. Begin Integration
  #########################################################################################
  .logDiffTime("Prepping Interation Data", tstart, verbose = verbose, logFile = logFile)
  
  #Clean Project For Parallel
  subProj <- ArchRProj
  subProj@imputeWeights <- SimpleList()
  
  #Gene Score Info
  geneDF <- .getFeatureDF(getArrowFiles(subProj), useMatrix)
  geneDF <- geneDF[geneDF$name %in% rownames(seuratRNA), , drop = FALSE]
  
  #Re-Index RNA
  splitGeneDF <- S4Vectors::split(geneDF, geneDF$seqnames)
  featureDF <- lapply(splitGeneDF, function(x){
    x$idx <- seq_len(nrow(x))
    return(x)
  }) %>% Reduce("rbind", .)
  dfParams <- data.frame(
    reduction = reduction
  )
  allChr <- unique(featureDF$seqnames)
  
  #Temp File Prefix
  tmpFile <- .tempfile()
  o <- suppressWarnings(file.remove(paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")))
  
  if(threads > 1){
    h5disableFileLocking()
  }
  
  rD <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  
  #Create Output Directory
  outDir1 <- getOutputDirectory(ArchRProj)
  outDir2 <- file.path(outDir1, "RNAIntegration")
  outDir3 <- file.path(outDir2, matrixName)
  dir.create(outDir1, showWarnings = FALSE)
  dir.create(outDir2, showWarnings = FALSE)
  dir.create(outDir3, showWarnings = FALSE)
  prevFiles <- list.files(outDir3, full.names = TRUE)
  prevFiles <- .suppressAll(file.remove(prevFiles))
  
  tstart <- Sys.time()
  
  threads2 <- max(ceiling(threads * 0.75), 1) #A Little Less here for now
  
  .logDiffTime(paste0("Computing Integration in ", length(blockList), " Integration Blocks!"), tstart, verbose = verbose, logFile = logFile)
  
  #Integration
  dfAll <- .safelapply(seq_along(blockList), function(i){
    
    prefix <- sprintf("Block (%s of %s) :", i , length(blockList))
    
    .logDiffTime(sprintf("%s Computing Integration", prefix), tstart, verbose = verbose, logFile = logFile)
    blocki <- blockList[[i]]
    
    #Subset ATAC
    subProj@cellColData <- subProj@cellColData[blocki$ATAC, ]
    subProj@sampleColData <- subProj@sampleColData[unique(subProj$Sample),,drop=FALSE]
    
    #Subset RNA
    subRNA <- seuratRNA[, blocki$RNA]
    
    #Subet RNA
    subRNA <- subRNA[rownames(subRNA) %in% geneDF$name, ]
    
    ##############################################################################################
    #1. Create Seurat RNA and Normalize
    ##############################################################################################
    .logDiffTime(sprintf("%s Identifying Variable Genes", prefix), tstart, verbose = verbose, logFile = logFile)
    subRNA <- FindVariableFeatures(object = subRNA, nfeatures = nGenes, verbose = FALSE)
    subRNA <- ScaleData(object = subRNA, verbose = FALSE)
    if(is.null(genesUse)){
      genesUse <- VariableFeatures(object = subRNA)
    }
    
    ##############################################################################################
    #2. Get Gene Score Matrix and Create Seurat ATAC
    ##############################################################################################
    .logDiffTime(sprintf("%s Getting GeneScoreMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
    mat <- .getPartialMatrix(
      getArrowFiles(subProj), 
      featureDF = geneDF[geneDF$name %in% genesUse,], 
      threads = 1,
      cellNames = subProj$cellNames,
      useMatrix = useMatrix,
      verbose = FALSE
    )
    rownames(mat) <- geneDF[geneDF$name %in% genesUse, "name"]
    .logThis(mat, paste0("GeneScoreMat-Block-",i), logFile=logFile)
    
    #Impute Matrix (its already scaled internally in ArrowFiles)
    if(useImputation){
      .logDiffTime(sprintf("%s Imputing GeneScoreMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
      imputeParams <- list()
      imputeParams$ArchRProj <- subProj
      imputeParams$randomSuffix <- TRUE
      imputeParams$reducedDims <- reducedDims
      imputeParams$dimsToUse <- dimsToUse
      imputeParams$scaleDims <- scaleDims
      imputeParams$corCutOff <- corCutOff
      imputeParams$threads <- 1
      imputeParams$logFile <- logFile
      subProj <- suppressMessages(do.call(addImputeWeights, imputeParams))
      mat <- suppressMessages(imputeMatrix(mat = mat, imputeWeights = getImputeWeights(subProj), verbose = FALSE, logFile = logFile))
      o <- suppressWarnings(file.remove(unlist(getImputeWeights(subProj)[[1]]))) #Clean Up Space
      .logThis(mat, paste0("GeneScoreMat-Block-Impute-",i), logFile=logFile)
    }
    
    #Log-Normalize 
    mat <- log(mat + 1) #use natural log
    seuratATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    
    #Clean Memory
    rm(mat)
    
    #Set Default Assay
    DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = FALSE)
    
    ##############################################################################################
    #3. Transfer Anchors  
    ############################################################################################## 
    .logDiffTime(sprintf("%s Seurat FindTransferAnchors", prefix), tstart, verbose = verbose, logFile = logFile)
    transferAnchors <- .retryCatch({ #This sometimes can crash in mclapply so we can just add a re-run parameter
      gc()
      Seurat::FindTransferAnchors(
        reference = subRNA, 
        query = seuratATAC, 
        reduction = reduction, 
        features = genesUse,
        verbose = FALSE,
        ...
      )
    }, maxAttempts = 2, logFile = logFile)
    .logThis(paste0(utils::capture.output(transferAnchors),collapse="\n"), paste0("transferAnchors-",i), logFile=logFile)
    
    ##############################################################################################
    #4. Transfer Data
    ##############################################################################################
    rDSub <- rD[colnames(seuratATAC),,drop=FALSE]
    .logThis(rDSub, paste0("rDSub-", i), logFile = logFile)
    transferParams$anchorset <- transferAnchors
    transferParams$weight.reduction <- CreateDimReducObject(
      embeddings = rDSub, 
      key = "LSI_", 
      assay = DefaultAssay(seuratATAC)
    )
    transferParams$verbose <- FALSE
    transferParams$dims <- seq_len(ncol(rDSub))
    transferParams$k.weight <- k.weight
    
    #Group
    .logDiffTime(sprintf("%s Seurat TransferData Cell Group Labels", prefix), tstart, verbose = verbose, logFile = logFile)
    transferParams$refdata <- subRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)
    
    #RNA Names
    .logDiffTime(sprintf("%s Seurat TransferData Cell Names Labels", prefix), tstart, verbose = verbose, logFile = logFile)
    transferParams$refdata <- colnames(subRNA)
    rnaLabels2 <- do.call(Seurat::TransferData, transferParams)[,1]
    
    if(addToArrow){
      .logDiffTime(sprintf("%s Seurat TransferData GeneMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
      transferParams$refdata <- GetAssayData(subRNA, assay = "RNA", slot = "data")
      gc()
      matchedRNA <- do.call(Seurat::TransferData, transferParams)
      matchedRNA <- matchedRNA@data
    }
    
    #Match results
    matchDF <- DataFrame(
      cellNames = colnames(seuratATAC), 
      predictionScore = rnaLabels$prediction.score.max,
      predictedGroup = rnaLabels$predicted.id,
      predictedCell = rnaLabels2
    )
    rownames(matchDF) <- matchDF$cellNames
    
    .logDiffTime(sprintf("%s Saving TransferAnchors Joint CCA", prefix), tstart, verbose = verbose, logFile = logFile)
    jointCCA <- DataFrame(transferAnchors@object.list[[1]]@reductions$cca@cell.embeddings)
    jointCCA$Assay <- ifelse(endsWith(rownames(jointCCA), "_reference"), "RNA", "ATAC")
    jointCCA$Group <- NA
    jointCCA$Score <- NA
    jointCCA[paste0(colnames(subRNA), "_reference"), "Group"] <- subRNA$Group
    jointCCA[paste0(matchDF$cellNames, "_query"), "Group"] <- matchDF$predictedGroup
    jointCCA[paste0(matchDF$cellNames, "_query"), "Score"] <- matchDF$predictionScore
    .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
    
    #Clean Memory
    rm(transferParams, transferAnchors)
    gc()
    
    ##############################################################################################
    #5. Add To Temp Hdf5
    ##############################################################################################
    
    if(addToArrow){
      
      .logDiffTime(sprintf("%s Transferring Paired RNA to Temp File", prefix), tstart, verbose = verbose, logFile = logFile)
      
      #Quickly Write to A Temp Hdf5 File Split By Sample to Then Enable Writing to Each Arrow File
      
      tmpFilei <- paste0(tmpFile, "-IntegrationBlock-", i, ".h5")
      o <- h5createFile(tmpFilei)
      sampleNames <- getCellColData(subProj, "Sample")[matchDF$cellNames, ]
      uniqueSamples <- unique(sampleNames)
      matchedRNA <- .safeSubset( #If Rownames disappeared this will catch that!
        mat = matchedRNA, 
        subsetRows = paste0(featureDF$name), 
        subsetCols = matchDF$cellNames
      )
      
      for(z in seq_along(uniqueSamples)){
        
        mat <- matchedRNA[, which(sampleNames == uniqueSamples[z]), drop = FALSE]
        Group <- uniqueSamples[z]
        
        o <- tryCatch({h5delete(tmpFilei, paste0(Group))}, error = function(x){})
        o <- h5createGroup(tmpFilei, paste0(Group))
        
        #Convert Columns to Rle
        j <- Rle(findInterval(seq(mat@x)-1, mat@p[-1]) + 1)
        
        #Info
        lengthRle <- length(j@lengths)
        lengthI <- length(mat@i)
        
        #Create Data Set
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/i"), storage.mode = "integer", 
                                          dims = c(lengthI, 1), level = 0))
        
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jLengths"), storage.mode = "integer", 
                                          dims = c(lengthRle, 1), level = 0))
        
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jValues"), storage.mode = "integer", 
                                          dims = c(lengthRle, 1), level = 0))
        
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group, "/x"), storage.mode = "double", 
                                          dims = c(lengthI, 1), level = 0))
        
        #Write Data Set
        o <- .suppressAll(h5write(obj = mat@i + 1, file = tmpFilei, name = paste0(Group,"/i")))
        o <- .suppressAll(h5write(obj = j@lengths, file = tmpFilei, name = paste0(Group,"/jLengths")))
        o <- .suppressAll(h5write(obj = j@values, file = tmpFilei, name = paste0(Group,"/jValues")))
        o <- .suppressAll(h5write(obj = mat@x, file = tmpFilei, name = paste0(Group, "/x")))
        o <- .suppressAll(h5write(obj = colnames(mat), file = tmpFilei, name = paste0(Group, "/cellNames")))
        #Row Names is always the same
        
      }
      
      rm(matchedRNA, mat, j)
      
    }
    
    .logDiffTime(sprintf("%s Completed Integration", prefix), tstart, verbose = verbose, logFile = logFile)
    
    gc()
    
    matchDF$Block <- Rle(i)
    matchDF
    
  }, threads = threads2) %>% Reduce("rbind", .)
  
  ##############################################################################################
  #5. Plot UMAPs for Co-Embeddings from CCA
  ##############################################################################################
  if(plotUMAP){
    
    for(i in seq_along(blockList)){
      
      o <- tryCatch({
        
        prefix <- sprintf("Block (%s of %s) :", i , length(blockList))
        
        .logDiffTime(sprintf("%s Plotting Joint UMAP", prefix), tstart, verbose = verbose, logFile = logFile)
        
        jointCCA <- readRDS(file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
        
        set.seed(1) # Always do this prior to UMAP
        UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors = 40, min_dist = 0.4, metric="cosine", verbose=FALSE))
        UMAPParams$X <- as.data.frame(jointCCA[, grep("CC_", colnames(jointCCA))])
        UMAPParams$ret_nn <- FALSE
        UMAPParams$ret_model <- FALSE
        UMAPParams$n_threads <- 1
        uwotUmap <- tryCatch({
          do.call(uwot::umap, UMAPParams)
        }, error = function(e){
          errorList <- UMAPParams
          .logError(e, fn = "uwot::umap", info = prefix, errorList = errorList, logFile = logFile)
        })
        
        #Add UMAP and Save Again
        jointCCA$UMAP1 <- uwotUmap[,1]
        jointCCA$UMAP2 <- uwotUmap[,2]
        .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
        
        p1 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = jointCCA$Assay,
          randomize = TRUE, 
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        p2 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = jointCCA$Group, 
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by scRNA Group"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        pdf(file.path(outDir3, paste0("Save-Block", i,"-JointCCA-UMAP.pdf")), width = 12, height = 6, useDingbats = FALSE)
        ggAlignPlots(p1,p2,type="h")
        dev.off()
        
      }, error = function(e){
        
      })
      
    }
    
  }
  
  ##############################################################################################
  #6. Read sub-matrices and store in ArrowFiles
  ##############################################################################################
  
  if(addToArrow){
    
    .logDiffTime("Transferring Data to ArrowFiles", tstart, verbose = verbose, logFile = logFile)
    
    matrixName <- .isProtectedArray(matrixName)
    
    integrationFiles <- paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")
    
    if(!all(file.exists(integrationFiles))){
      .logMessage("Something went wrong with integration as not all temporary files containing integrated RNA exist!", logFile = logFile)
      stop("Something went wrong with integration as not all temporary files containing integrated RNA exist!")
    }
    
    h5list <- .safelapply(seq_along(integrationFiles), function(x){
      h5ls(integrationFiles[x])
    }, threads = threads)
    
    ArrowFiles <- getArrowFiles(ArchRProj)
    allSamples <- names(ArrowFiles)
    
    o <- .safelapply(seq_along(allSamples), function(y){
      
      sample <- allSamples[y]
      
      prefix <- sprintf("%s (%s of %s)", sample, y, length(ArrowFiles))
      
      .logDiffTime(sprintf("%s Getting GeneIntegrationMatrix From TempFiles!", prefix), tstart, verbose = verbose, logFile = logFile)
      
      sampleIF <- lapply(seq_along(h5list), function(x){
        if(any(h5list[[x]]$group==paste0("/",sample))){
          integrationFiles[x]
        }else{
          NULL
        }
      }) %>% unlist
      
      sampleMat <- lapply(seq_along(sampleIF), function(x){
        
        cellNames <- .h5read(sampleIF[x], paste0(sample, "/cellNames"))
        
        mat <- sparseMatrix(
          i = .h5read(sampleIF[x], paste0(sample, "/i"))[,1], 
          j = as.vector(
            Rle(
              .h5read(sampleIF[x], paste0(sample, "/jValues"))[,1], 
              .h5read(sampleIF[x], paste0(sample, "/jLengths"))[,1]
            )
          ), 
          x = .h5read(sampleIF[x], paste0(sample, "/x"))[,1],
          dims = c(nrow(featureDF), length(cellNames))
        )
        colnames(mat) <- cellNames
        
        mat
        
      }) %>% Reduce("cbind", .)
      
      sampleMat@x <- exp(sampleMat@x) - 1 #Back To Counts
      sampleMat <- .normalizeCols(sampleMat, scaleTo = scaleTo) #Scale to 10,000
      sampleMat <- drop0(sampleMat) # Drop 0's
      rownames(sampleMat) <- paste0(featureDF$name)
      sampleMat <- sampleMat[,ArchRProj$cellNames[BiocGenerics::which(ArchRProj$Sample == sample)], drop = FALSE]
      
      ######################################
      # Initialize SP Mat Group
      ######################################
      o <- .createArrowGroup(ArrowFile = ArrowFiles[sample], group = matrixName, force = force)
      
      o <- .initializeMat(
        ArrowFile = ArrowFiles[sample],
        Group = matrixName,
        Class = "double",
        Units = "NormCounts",
        cellNames = colnames(sampleMat),
        params = dfParams,
        featureDF = featureDF,
        force = force
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictionScore"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictionScore")
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedGroup"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictedGroup")
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedCell"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictedCell")
      )
      
      .logDiffTime(sprintf("%s Adding GeneIntegrationMatrix to ArrowFile!", prefix), tstart, verbose = verbose, logFile = logFile)
      
      for(z in seq_along(allChr)){
        
        chrz <- allChr[z]
        
        .logDiffTime(sprintf("Adding GeneIntegrationMatrix to %s for Chr (%s of %s)!", sample, z, length(allChr)), tstart, verbose = FALSE, logFile = logFile)
        
        idz <- BiocGenerics::which(featureDF$seqnames %bcin% chrz)
        matz <- sampleMat[idz, ,drop=FALSE]
        stopifnot(identical(paste0(featureDF$name[idz]), paste0(rownames(matz))))
        
        #Write sparseMatrix to Arrow File!
        o <- .addMatToArrow(
          mat = matz, 
          ArrowFile = ArrowFiles[sample], 
          Group = paste0(matrixName, "/", chrz), 
          binarize = FALSE,
          addColSums = TRUE,
          addRowSums = TRUE,
          addRowVarsLog2 = TRUE,
          logFile = logFile
        )
        
        #Clean Memory
        rm(matz)
        
        if(z %% 3 == 0 | z == length(allChr)){
          gc()
        }
        
      }
      
      0
      
    }, threads = threads)
    
    o <- suppressWarnings(file.remove(integrationFiles))
    
  }
  
  .logDiffTime("Completed Integration with RNA Matrix", tstart, verbose = verbose, logFile = logFile)
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictedCell,
    name = nameCell,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictedGroup,
    name = nameGroup,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictionScore,
    name = nameScore,
    force = TRUE
  )
  
  .endLogging(logFile = logFile)
  
  return(ArchRProj)
  
}

environment(my_addGeneIntegrationMatrix) <- asNamespace('ArchR')

transform_p2g <- function(p2g) {
  p2g$idxATAC <- (metadata(p2g)$peakSet %>% as.character)[p2g$idxATAC]
  p2g$idxRNA <- metadata(p2g)$geneSet$name[p2g$idxRNA]
  
  return(p2g)
}

my_loopTracks <- function(
  loops = NULL, 
  region = NULL, 
  title = "LoopTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = 9,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
){
  
  getArchDF <- function(lp, r = 100){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)
    return(df)
  }
  
  if(!is.null(loops)){
    
    if(is(loops, "GRanges")){
      loops <- SimpleList(Loops = loops)
    }else if(ArchR:::.isGRList(loops)){
    }else{
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }
    
    valueMin <- min(unlist(lapply(loops, function(x) min(x$value))))
    valueMax <- max(unlist(lapply(loops, function(x) max(x$value))))
    
    loopO <- lapply(seq_along(loops), function(x){
      subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within") 
      if(length(subLoops)>0){
        dfx <- getArchDF(subLoops)
        dfx$name <- Rle(paste0(names(loops)[x]))
        dfx
      }else{
        NULL
      }
    }) %>% Reduce("rbind",.)
    loopO$name <- factor(x = loopO$name, levels = names(loops))
    .logThis(loopO, "loopO", logFile = logFile)
    
    testDim <- tryCatch({
      if(is.null(loopO)){
        FALSE
      }
      if(nrow(loopO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })
    
    if(testDim){
      
      loopO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
      }
      
      p <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line() +
        facet_grid(name ~ ., drop = F) +
        ylab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax)) +
        theme(legend.text = element_text(size = baseSize)) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
              legend.box.background = element_rect(color = NA)) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3))
      
    }else{
      
      #create empty plot
      df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
      p <- ggplot(data = df, aes(start, end)) + 
        geom_point() +
        facet_grid(facet~.) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
      
    }
    
  }else{
    
    #create empty plot
    df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
    
  }
  
  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  
  if(!is.ggplot(p)){
    .logError("loopTracks is not a ggplot!", fn = ".loopTracks", info = "", errorList = NULL, logFile = logFile)
  }
  
  return(p)
  
}

environment(my_loopTracks) <- asNamespace('ArchR')


# genes, regions are mutex
# if select use_matrix="GeneIntegrationMatrix", rna_files must be provided, or Arch_proj has GeneIntegrationMatrix
plot_peak2gene_link <- function(Arch_proj, genes = NULL, regions = NULL, upstream = 1e5, downstream = 1e5, 
                                track_sizes = 3, title_size = 16, gene_size = 6, 
                                outdir = "peak2gene", groupby = "tissue", use_groups = NULL, use_matrix = "GeneIntegrationMatrix", 
                                rna_sobjs = NULL, cor_cut_off = 0.6, FDR_cut_off = 1e-4, promoter_length = c(2000, 100), promoter_only = T, 
                                save_p2g_mat = F, group_rna = "group", logFile = createLogFile("plot_peak2gene_link")) {
  ArchR:::.startLogging(logFile = logFile)
  
  # construction of peak2gene info
  if (is.null(use_groups)) {
    use_groups <- getCellColData(Arch_proj, select = groupby, drop = T) %>% unique
  } else {
    use_groups <- use_groups[use_groups %in% getCellColData(Arch_proj, select = groupby, drop = T) %>% unique]
  }
  
  cat("analyzing groups:", use_groups, "\n")
  
  dir.create(outdir, showWarnings = F, recursive = T)
  loops_by_groups <- lapply(use_groups, function(group) {
    cat(group, "\n")
    proj_one_group_path <- paste( outdir, gsub(" ", "_", group), sep = "/")
    if( dir.exists(proj_one_group_path) ){
      proj_one_group <- loadArchRProject(proj_one_group_path, showLogo = F)
    } else {
      proj_one_group <- Arch_proj[getCellColData(Arch_proj, groupby, T) == group]
      
      # first save and then perform downstream analysis, this will cause changes in getSampleColData(proj_sub, )[, "ArrowFiles"]
      # however, while computing LSI, ArchR first subsets peak/tile to highly variable or top ones, which is according to arrowfiles (not subsetted)
      # so this will cause some problems, as drop cells should be implemented
      # error code as follows: totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
      proj_one_group <- saveArchRProject(proj_one_group, proj_one_group_path, dropCells = T)
      
      proj_one_group <- addIterativeLSI(proj_one_group, useMatrix = "PeakMatrix", 
                                        iterations = 1, name = "peakLSI", varFeatures = 50000, force = T, verbose = F)
      ArchR:::.logThis(x = getReducedDims(proj_one_group, reducedDims = "peakLSI"), name = "peakLSI", logFile = logFile)
      
      proj_one_group <- addUMAP(ArchRProj = proj_one_group, reducedDims = "peakLSI", name = "peakUMAP", 
                                nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T, verbose = F)
      ArchR:::.logThis(x = getEmbedding(proj_one_group, embedding = "peakUMAP"), name = "peakUMAP", logFile = logFile)
      # p1 <- plotEmbedding(proj_one_group, embedding = "peakUMAP", colorBy = "cellColData", name = groupby, size = 1)
      
      # determine used matrix for peak2gene correlation calculation
      proj_one_group <- addGeneScoreMatrix(proj_one_group, force = T)
      tmp <- getMatrixFromProject(proj_one_group, useMatrix = "GeneScoreMatrix") %>% assay
      ArchR:::.logThis(x = tmp, name = "GeneScoreMatrix", logFile = logFile)
      ArchR:::.logThis(x = tmp@x, name = "GeneScoreMatrix_x", logFile = logFile)
      
      if(use_matrix == "GeneIntegrationMatrix") {
        if(is.null(rna_sobjs)) {
          assert_that("GeneIntegrationMatrix" %in% getAvailableMatrices(proj_one_group), msg = "GeneIntegrationMatrix 
                      must exist for ArchR project, if rna seurat objects are not provided")
        } else {
          sobj_group <- rna_sobjs[[group]]
          cat( sprintf("Integration by %d RNA cells and %d ATAC cells\n", dim(sobj_group)[2], length( getCellNames(proj_one_group) ) ) )
          
          proj_one_group <- my_addGeneIntegrationMatrix(proj_one_group, useMatrix = "GeneScoreMatrix",
                                                        reducedDims = "peakLSI", seRNA = sobj_group,
                                                        dimsToUse = 1:(min( 0.4 * length( getCellNames(proj_one_group) ), 30 ) %>% floor),
                                                        addToArrow = TRUE,
                                                        force= TRUE,
                                                        groupRNA = group_rna, # TODO 这里有没有更好的写法
                                                        nameCell = "predictedCell",
                                                        nameGroup = "predictedGroup",
                                                        nameScore = "predictedScore", 
                                                        plotUMAP = F, 
                                                        useImputation = F, 
                                                        npcs = min( 0.4 * ncol(sobj_group), 0.4 * length( getCellNames(proj_one_group) ), 30 ) %>% floor, 
                                                        dims = 1:(min( 0.4 * ncol(sobj_group), 0.4 * length( getCellNames(proj_one_group) ), 30 ) %>% floor),
                                                        k.filter = min( 0.4 * ncol(sobj_group), 0.4 * length( getCellNames(proj_one_group) ), 200 ) %>% floor, 
                                                        k.score = min( 0.4 * ncol(sobj_group), 0.4 * length( getCellNames(proj_one_group) ), 30 ) %>% floor, 
                                                        k.weight = min( 0.4 * ncol(sobj_group), 0.4 * length( getCellNames(proj_one_group) ), 30 ) %>% floor)
          tmp <- getMatrixFromProject(proj_one_group, useMatrix = "GeneIntegrationMatrix") %>% assay
          ArchR:::.logThis(x = tmp, name = "GeneIntegrationMatrix", logFile = logFile)
          ArchR:::.logThis(x = tmp@x, name = "GeneIntegrationMatrix_x", logFile = logFile)
        }
      }
      
      proj_one_group <- addPeak2GeneLinks(proj_one_group, 
                                          useMatrix = use_matrix, 
                                          reducedDims = "peakLSI", 
                                          k = min( 0.4 * length( getCellNames(proj_one_group) ), 100 ) %>% floor)
      
      # change saving format for group RNA and ATAC data
      seRNA <- readRDS(paste(proj_one_group_path, "Peak2GeneLinks/seRNA-Group-KNN.rds", sep = "/"))
      seATAC <- readRDS(paste(proj_one_group_path, "Peak2GeneLinks/seATAC-Group-KNN.rds", sep = "/"))
      
      elementMetadata(seATAC)[["name"]] <- rowRanges(seATAC) %>% as.character
      
      colnames(seRNA) <- paste0("pseudo_cell_", seq_len( ncol(seRNA) ) )
      colnames(seATAC) <- paste0("pseudo_cell_", seq_len( ncol(seATAC) ) )
      
      sRNA <- summarizedExperiment2Seurat(seRNA, assay = "RawRNA", do_norm = T, rename_assay = "RNA")
      sATAC <- summarizedExperiment2Seurat(seATAC, assay = "RawATAC", do_norm = T, rename_assay = "ATAC")
      
      saveRDS(sRNA, file = paste(proj_one_group_path, "Peak2GeneLinks/seuratRNA-Group-KNN.rds", sep = "/") )
      saveRDS(sATAC, file = paste(proj_one_group_path, "Peak2GeneLinks/seuratATAC-Group-KNN.rds", sep = "/") )
      
      saveArchRProject(proj_one_group, proj_one_group_path)
      
      p2g <- getPeakSet(proj_one_group) %>% metadata %>% .$Peak2GeneLinks
      p2g <- p2g[!is.na(p2g$Correlation), ]
      p2g <- transform_p2g(p2g)
      saveRDS(p2g, paste( proj_one_group_path, "Peak2GeneLinks/peak2gene.rds", sep = "/") )
      
      p2g <- as.data.frame(p2g)
      colnames(p2g)[c(1, 2)] <- c("peak", "gene")
      
      write.table(p2g, paste(proj_one_group_path, "Peak2GeneLinks/peak2gene.txt", sep = "/"), sep='\t',
                  quote = F, row.names = F, col.names = T)
    }
    
    loops <- getPeak2GeneLinks(proj_one_group, 
                               corCutOff = cor_cut_off, 
                               FDRCutOff = FDR_cut_off)
    ArchR:::.logThis(x = loops, name = "loops", logFile = logFile)
    
    if(save_p2g_mat) {
      p2g_heatmap_mat <- plotPeak2GeneHeatmap(proj_one_group, corCutOff = cor_cut_off, FDRCutOff = FDR_cut_off, 
                                              k = min( length(loops[[1]]) * 0.2, 25) %>% floor, groupBy = groupby, 
                                              returnMatrices = T, verbose = F)
      saveRDS(p2g_heatmap_mat, paste(proj_one_group_path, "Peak2GeneLinks/p2g_heatmap_mat.rds", sep = "/") )
      
      write.table(p2g_heatmap_mat$Peak2GeneLinks, file = paste(proj_one_group_path, "Peak2GeneLinks/filtered_peak2gene.txt", sep = "/"), 
                  sep = "\t", quote = F, col.names = NA)
    }
    
    return(loops[[1]])
  })
  
  names(loops_by_groups) <- use_groups
  
  # plotting
  
  geneAnnotation <- getGeneAnnotation(Arch_proj)
  
  if (!is.null(genes)) {
    regions <- geneAnnotation$genes
    # only keep genes with gene annotation
    genes <- genes[tolower(genes) %in% tolower(mcols(regions)$symbol)]
    
    regions <- regions[which(tolower(mcols(regions)$symbol) %in% tolower(genes))]
    regions <- regions[order(match(tolower(mcols(regions)$symbol), tolower(genes)))]
    # print(regions)
    regions <- resize(regions, 1, "start")
    # strand(region) <- "*"
    if(promoter_only) {
      promoter_regions <- extendGR(regions, upstream = promoter_length[1], downstream = promoter_length[2])
    }
    regions <- extendGR(regions, upstream = upstream, downstream = downstream)
    names(promoter_regions) <- genes
    names(regions) <- genes
  }
  assert_that(!is.null(regions), msg = "regions or genes must be provided")
  
  trackfig_regions <- lapply(seq_along(regions), function(i) {
    # subset loops according to region
    # Note: r has rather weird variable scope, but is exactly what we need here (function don't change global variables)
    if(promoter_only & !is.null(genes)) {
      loops_by_groups <- lapply(loops_by_groups, function(loops) {
        promoter_region <- promoter_regions[names(regions)[i]]
        
        loops <- loops[(start(loops) > start(promoter_region) & start(loops) < end(promoter_region)) |
                         (end(loops) > start(promoter_region) & end(loops) < end(promoter_region))]
        # print(loops)
        loops
      })
    }
    loop_track <- my_loopTracks(loops = loops_by_groups, 
                                region = regions[i], 
                                facetbaseSize = title_size,
                                hideX = TRUE,
                                hideY = TRUE,
                                title = "Loops",
                                logFile = NULL) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    
    # add gene track
    gene_track <- ArchR:::.geneTracks(geneAnnotation = geneAnnotation, region = regions[i], facetbaseSize = title_size, title = "Genes", 
                                      labelSize = gene_size, logFile = createLogFile(name = "plot_peak2gene_links")) + 
      theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    
    p_list <- SimpleList(loop_track, gene_track)
    # TODO viewpoint
    
    # merge all tracks
    ArchR:::ggAlignPlots(plotList = p_list, sizes=c(3*length(use_groups), 3), draw = FALSE)
    
  })
  
  if (!is.null(genes))
    names(trackfig_regions) <- genes
  
  ArchR:::.endLogging(logFile = logFile)
  trackfig_regions
}
