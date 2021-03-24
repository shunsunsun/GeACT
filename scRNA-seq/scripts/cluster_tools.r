
# set lib for UMAP
if(! grepl("/home/tianf/lustre/Tools/Python-3.6.4/install/lib", Sys.getenv("LD_LIBRARY_PATH"))) {
  #print("Set LD_LIBRARY_PATH for Python lib.")
  Sys.setenv(LD_LIBRARY_PATH=paste0("/home/tianf/lustre/Tools/Python-3.6.4/install/lib:",Sys.getenv("LD_LIBRARY_PATH")))
}

do_addMeta <- function(expr) {
  batch <- gsub("_[0-9]+$", "", colnames(expr@raw.data))
  names(batch) <- colnames(expr@raw.data)
  bigBatch <- gsub("-.*", "", colnames(expr@raw.data))
  names(bigBatch) <- colnames(expr@raw.data)
  outer_id <- ceiling(as.numeric(gsub(".*_", "", colnames(expr@raw.data))) / 8); outer_id <- factor(outer_id)
  names(outer_id) <- colnames(expr@raw.data)
  inner_id <- as.numeric(gsub(".*_", "", colnames(expr@raw.data))) %% 8; inner_id[inner_id==0] <- 8; inner_id <- factor(inner_id, levels = 1:8)
  names(inner_id) <- colnames(expr@raw.data)
  bigInner_id <- ifelse(as.numeric(inner_id)<=4, "1-4", "5-8")
  names(bigInner_id) <- colnames(expr@raw.data)
  
  expr <- AddMetaData(object = expr, metadata = percent.mito, col.name = "percent.mito")
  expr <- AddMetaData(object = expr, metadata = batch, col.name = "batch")
  expr <- AddMetaData(object = expr, metadata = bigBatch, col.name = "bigBatch")
  #expr <- AddMetaData(object = expr, metadata = people, col.name = "people")
  expr <- AddMetaData(object = expr, metadata = outer_id, col.name = "outer_id")
  expr <- AddMetaData(object = expr, metadata = inner_id, col.name = "inner_id")
  expr <- AddMetaData(object = expr, metadata = bigInner_id, col.name = "bigInner_id")
  
  return(expr)
}


do_batchDebug <- function() {
  p <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + geom_point() + theme_bw() + 
    theme(legend.background = element_blank(), legend.box.background = element_rect()) #+ 
  #theme(panel.grid = element_blank(), legend.position = c(1.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=2,byrow=TRUE))
  print(p)
  
  p <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = batch)) + geom_point(alpha = 0.8, show.legend = F, size = 0.6) + theme_bw() + 
    theme(legend.background = element_blank(), legend.box.background = element_rect(fill=alpha('white', 0.4))) #+ 
  #theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
  print(p)
  
  ###
  for(i in levels(cellMetaData$bigBatch)) {
    tmp <- cellMetaData[cellMetaData$bigBatch==i, c("batch", "cluster")]
    tmp$batch <- factor(tmp$batch, levels = unique(tmp$batch))
    tmp <- as.data.frame(table(tmp))
    #p <- ggplot(tmp, aes(x = cluster, y = Freq, fill = batch)) + geom_bar(stat = "identity", show.legend = T)
    p <- ggplot(tmp, aes(x = cluster, y = Freq, fill = batch)) + geom_bar(stat = "identity", position="fill", show.legend = T) + 
      scale_y_continuous(expand = c(0, 0))
    print(p)
    p <- ggplot(subset(cellMetaData, bigBatch==i), aes(x = tSNE_1, y = tSNE_2, color = batch)) + geom_point(show.legend = T) + theme_bw() + 
      theme(legend.background = element_blank(), legend.box.background = element_rect(fill=alpha('white', 0.4))) #+ 
    #theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + guides(color = guide_legend(ncol=1,byrow=TRUE))
    print(p)
  }
  ###
  
  p <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = bigBatch)) + geom_point() + theme_bw() + 
    theme(legend.background = element_blank(), legend.box.background = element_rect()) #+ 
  print(p)
  p <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = outer_id)) + geom_point() + theme_bw() + 
    theme(legend.background = element_blank(), legend.box.background = element_rect(fill = alpha('white', 0.4))) + 
    theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + 
    guides(color = guide_legend(ncol=2,byrow=TRUE))
  print(p)
  p <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = inner_id)) + geom_point() + theme_bw() + 
    theme(legend.background = element_blank(), legend.box.background = element_rect(fill = alpha('white', 0.4))) + 
    theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + 
    guides(color = guide_legend(ncol=1,byrow=TRUE))
  print(p)
  p <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = bigInner_id)) + geom_point() + theme_bw() + 
    theme(legend.background = element_blank(), legend.box.background = element_rect(fill = NA)) + 
    theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + 
    guides(color = guide_legend(ncol=1,byrow=TRUE))
  print(p)
  
  # other features
  for(i in colnames(cellMetaData)[c(2,3,15:20)]) {
    gp <- ggplot(cellMetaData, aes(x = tSNE_1, y = tSNE_2, color = get(i))) + geom_point() + theme_bw() + 
      theme(legend.background = element_blank(), legend.box.background = element_rect(fill = alpha('white', 0.4))) + 
      theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + 
      scale_color_continuous(name = i)
    print(gp)
  }
  
  for(i in colnames(cellMetaData)[c(2,3,15:20)]) {
    gp <- ggplot(cellMetaData, aes(x = batch, y = get(i), fill = batch)) + geom_boxplot(show.legend = F, notch = T) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("Plate") + 
      theme(panel.grid = element_blank(), legend.position = c(0.01, 0.99), legend.justification = c("left", "top")) + 
      ylab(i)
    print(gp)
  }
  
  for(i in "ERCCratio") {
    gp <- ggplot(cellMetaData, aes(x = cluster, y = get(i), fill = cluster)) + geom_boxplot(show.legend = F, notch = T) + 
      ylab(i)
    print(gp)
  }
}


do_checkMarker <- function(x, marker, cell) {
  expr_sub <- x[marker, grep(cell, colnames(x))]
  rownames(expr_sub) <- marker
  #cat("Data dimension: ", dim(expr_sub), "\n")
  expr_ratio <- t(apply(expr_sub, 1, function(x) { y0 <- sum(x==0); y1 <- sum(x>0); y2 <- y1/(y0 + y1); return(c(y0, y1, y2)) }))
  colnames(expr_ratio) <- c("Non_expressed", "Expressed", "ratio")
  expr_ratio <- data.frame(bigBatch=cell, marker=rownames(expr_ratio), expr_ratio)
  return(expr_ratio)
}


do_clusterDebug <- function() {
  ## bigBatch level
  # outer
  for(i in unique(cellMetaData$bigBatch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
      geom_violin(aes(x = outer_id, y = mitoRatio), fill = "navy", alpha = 0.3) + 
      ggtitle(paste(i, "(outer barcode: mito)"))
    print(gp)
  }
  
  for(i in unique(cellMetaData$bigBatch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
      geom_bar(aes(x = outer_id, fill = cluster)) + 
      ggtitle(paste(i, "(outer barcode)")) + 
      scale_y_continuous(expand = c(0, 0)) + 
      ylab("Cell number")
    print(gp)
  }
  
  for(i in unique(cellMetaData$bigBatch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
      geom_bar(aes(x = outer_id, fill = cluster), position="fill") + 
      ggtitle(paste(i, "(outer barcode)")) + 
      scale_y_continuous(expand = c(0, 0)) + 
      ylab("Cell frequency")
    print(gp)
  }
  
  # inner
  for(i in unique(cellMetaData$bigBatch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
      geom_violin(aes(x = inner_id, y = mitoRatio), fill = "navy", alpha = 0.3) + 
      ggtitle(paste(i, "(inner barcode: mito)"))
    print(gp)
  }
  
  for(i in unique(cellMetaData$bigBatch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
      geom_bar(aes(x = inner_id, fill = cluster)) + 
      ggtitle(paste(i, "(inner barcode)")) + 
      scale_y_continuous(expand = c(0, 0)) + 
      ylab("Cell number")
    print(gp)
  }
  
  for(i in unique(cellMetaData$bigBatch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$bigBatch==i, ]) + 
      geom_bar(aes(x = inner_id, fill = cluster), position="fill") + 
      ggtitle(paste(i, "(inner barcode)")) + 
      scale_y_continuous(expand = c(0, 0)) + 
      ylab("Cell frequency")
    print(gp)
  }
  
  ## batch level
  for(i in unique(cellMetaData$batch)) {
    print(i)
    gp <- ggplot(cellMetaData[cellMetaData$batch==i, ]) + 
      geom_bar(aes(x = 1, y = sqrt(cleanReads), fill = mitoRatio, color = cluster), alpha = 1, stat = "identity", show.legend = T) + 
      coord_polar("x") + 
      geom_text(x=0, y=0, vjust = 2, aes(label = round(mitoRatio*100,2)), size=2) + 
      facet_grid(inner_id ~ outer_id, drop = F) + ggtitle(i) + 
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
      theme(strip.text = element_text(size = 6)) + 
      theme(plot.margin = unit(c(0.1,0.1,0,0), "cm")) + 
      scale_fill_gradientn(colours = c("white", "blue"))
    print(gp)
  }
}


do_literatureMarker <- function() {
  if(file.exists(paste0(OUT, "/mg_literature.txt"))) {
    print("Analyze the expression of genes in literature")
    mg_literature <- read.table(paste0(OUT, "/mg_literature.txt"), header = F, sep = "\t", stringsAsFactors = F)
    mg_literature_LS <- apply(mg_literature, 1, function(x) {
      ###x <- mg_literature[1, ]
      y0 <- strsplit(x[2], ", ")[[1]]
      y <- data.frame(type = x[1], gene = y0, stringsAsFactors = F, row.names = NULL)
    })
    mg_literature_DF <- do.call("rbind", mg_literature_LS)
    mg_literature_DF[! mg_literature_DF$gene %in% geneID$symbol, ]
    p <- DoHeatmap_new(object = expr, genes.use = mg_literature_DF$gene, genes.group = mg_literature_DF$type, slim.col.label = TRUE, remove.key = TRUE, cex.row = 8, do.plot = F)
    print(p)
    mg_literature_sub <- intersect(mg_literature_DF$gene, rownames(expr_data))
    mg_literature_sub_LS <- split(mg_literature_sub, ceiling(seq_along(mg_literature_sub) / 4))
    for(i in seq_along(mg_literature_sub_LS)) {
      x1 <- mg_literature_sub_LS[[i]]
      if(length(x1) < 4) { x1 <- c(x1, rep(x1[1], 4 - length(x1))) }
      FeaturePlot(object = expr, features.plot = x1, cols.use = c("grey", "blue"), reduction.use = "tsne", pt.size = 2)
    }
  } else {
    mg_literature_DF <- NULL
  }
  return(mg_literature_DF)
}


do_findMarker <- function(expr, only.pos = T, ncpu = 2) {
  library("parallel")
  cl <- makeCluster(getOption("cl.cores", ncpu))
  clusterExport(cl, list("expr", "FindMarkers"))
  expr.markers_LS <- parLapply(cl, levels(expr@ident), function(x) {
    if(sum(expr@ident == x) >= 3) {
      y <- FindMarkers(object = expr, ident.1 = x, test.use = "roc", only.pos = only.pos, min.pct = 0.25)
      y$cluster <- x
      y$gene <- rownames(y)
    } else {
      y <- NULL
    }
    return(y)
  })
  stopCluster(cl)
  expr.markers <- do.call("rbind", expr.markers_LS)
  return(expr.markers)
}


do_findSpecMarker <- function(expr, cl.excluded = NULL, ncpu = 2) {
  library("parallel")
  suppressMessages(library("stringr"))
  cts <- str_sort(names(table(expr@ident)[table(expr@ident) >= 3]), numeric = T)
  cts <- setdiff(cts, cl.excluded)
  
  cl <- makeCluster(getOption("cl.cores", ncpu))
  clusterExport(cl, varlist = c("expr.markers_ftd", "cts", "FindMarkers", "expr"), envir = environment())
  # pairwise comparison
  mk_LS <- parLapply(cl, cts, function(ct1) {
    genes <- subset(expr.markers_ftd, cluster == ct1, select = "gene", drop = T)
    mk_DF <- head(data.frame(expr.markers_ftd[, 1:7], ct1 = NA, ct2 = NA, gene = NA), n = 0)
    if(length(genes) == 0) {
      print("No candidate genes.")
    } else {
      for(ct2 in cts) {
        cat(">", ct1, ct2, "\n")
        if(ct1 == ct2) { next }
        tmp <- FindMarkers(object = expr, test.use = "roc", only.pos = TRUE, min.pct = 0.25, 
                           ident.1 = ct1, ident.2 = ct2, genes.use = genes)
        tmp_ftd <- subset(tmp, power>=0.4 & avg_logFC>=log(1.5))
        if(nrow(tmp_ftd) > 0) {
          tmp_ftd <- data.frame(tmp_ftd, ct1, ct2, gene = rownames(tmp_ftd), stringsAsFactors = F, row.names = NULL)
          mk_DF <- rbind(mk_DF, tmp_ftd)
        }
      }
    }
    return(mk_DF)
  })
  stopCluster(cl); rm(cl)
  names(mk_LS) <- cts
  # extract spec genes
  expr.markers_spec_LS <- lapply(cts, function(ct) {
    spec_gene <- names(table(mk_LS[[ct]]$gene))[table(mk_LS[[ct]]$gene) == (length(cts) - 1)]
    mk_sub <- subset(expr.markers_ftd, cluster == ct)
    mk_sub$spec <- mk_sub$gene %in% spec_gene
    return(mk_sub)
  })
  expr.markers_spec <- do.call("rbind", expr.markers_spec_LS)
  rownames(expr.markers_spec) <- NULL
  expr.markers_spec$cluster <- factor(expr.markers_spec$cluster, levels = cts)
  expr.markers_spec <- subset(expr.markers_spec, spec)
  return(expr.markers_spec)
}


# do_GOenrich_old <- function(ncpu = 2) {
#   library("enrichR")
#   library("stringr")
#   library("parallel")
#   #dbs <- listEnrichrDbs()
#   expr.markers_ftd_LS <- split(x = expr.markers_ftd$gene, f = expr.markers_ftd$cluster)
#   expr.markers_ftd_LS <- expr.markers_ftd_LS[str_sort(names(expr.markers_ftd_LS), numeric = T)]
#   cl <- makeCluster(getOption("cl.cores", ncpu))
#   clusterExport(cl, c("expr.markers_ftd_LS", "enrichr"), envir = environment())
#   #enriched_LS <- parLapply(cl, seq_along(expr.markers_ftd_LS), function(i) {
#   enriched_LS <- lapply(seq_along(expr.markers_ftd_LS), function(i) {
#     cat(">", i, "\n")
#     x <- expr.markers_ftd_LS[[i]]
#     y <- enrichr(genes = x, databases = "GO_Biological_Process_2018")
#     y[[1]]$cluster <- names(expr.markers_ftd_LS)[i]
#     return(y)
#   })
#   stopCluster(cl); rm(cl)
#   return(enriched_LS)
# }


do_GOenrich <- function(expr.markers_ftd, ncpu = 1, do.filtering = T) {
  suppressMessages(library("stringr"))
  library("parallel")
  
  # read GO anno
  geneID2GO <- readMappings("~/lustre/06-Human_cell_atlas/Data/GO/goa_human_gene2GO.txt")
  geneNames <- names(geneID2GO)
  
  expr.markers_ftd_LS <- split(x = expr.markers_ftd$gene, f = expr.markers_ftd$cluster)
  expr.markers_ftd_LS <- expr.markers_ftd_LS[str_sort(names(expr.markers_ftd_LS), numeric = T)]
  
  cl <- makeCluster(getOption("cl.cores", ncpu))
  clusterExport(cl, c("expr.markers_ftd_LS"), envir = environment())
  
  enriched_LS <- parLapply(cl, seq_along(expr.markers_ftd_LS), function(i) {
    cat(">", i, "\n")
    suppressMessages(library("topGO"))
    geneInput <- expr.markers_ftd_LS[[i]]
    cat(geneInput, sep = "\n")
    geneList <- factor(as.integer(geneNames %in% geneInput))
    names(geneList) <- geneNames
    
    sampleGOdata <- new("topGOdata", description = "GOenrich", ontology = "BP",
                        allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    
    allGO <- usedGO(object = sampleGOdata)
    allRes <- GenTable(sampleGOdata, resultFisher, topNodes = length(allGO), numChar = 500)
    colnames(allRes)[6] <- "p_value"
    allRes[allRes$p_value=="< 1e-30", "p_value"] <- 1e-30
    allRes$p_value <- as.numeric(allRes$p_value)
    # FDR
    allRes$q_value <- p.adjust(p = allRes$p_value, method = "fdr")
    # filtering
    if(do.filtering) {
      allRes <- subset(allRes, p_value < 0.05)
    }
    # add gene list
    go2gene_LS <- genesInTerm(object = sampleGOdata, whichGO = allRes$GO.ID)
    go2gene_AR <- sapply(go2gene_LS, function(x) {
      x <- x[x %in% geneInput]
      y <- paste(x, collapse = ",")
    })
    allRes$Genes <- go2gene_AR
    allRes$cluster <- rep(names(expr.markers_ftd_LS)[i], nrow(allRes))
    return(allRes)
  })
  stopCluster(cl); rm(cl)
  return(enriched_LS)
}


do_regionMarker <- function(expr, cellMetaData, left1, right1, top1, bottom1, left2 = NULL, right2, top2, bottom2) {
  expr_test <- expr
  regionCL <- subset(cellMetaData, tSNE_1 > left1 & tSNE_1 < right1 & tSNE_2 > bottom1 & tSNE_2 < top1, "cell", drop = T)
  if(is.null(left2)) {
    expr_test@ident <- factor(expr_test@ident, levels = c(levels(expr_test@ident), 888))
    expr_test@ident[names(expr_test@ident) %in% regionCL] <- 888
    regionMk <- FindMarkers(object = expr_test, ident.1 = 888, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
  } else {
    regionCR <- subset(cellMetaData, tSNE_1 > left2 & tSNE_1 < right2 & tSNE_2 > bottom2 & tSNE_2 < top2, "cell", drop = T)
    expr_test@ident <- factor(expr_test@ident, levels = c(levels(expr_test@ident), 888, 999))
    expr_test@ident[names(expr_test@ident) %in% regionCL] <- 888
    expr_test@ident[names(expr_test@ident) %in% regionCR] <- 999
    regionMk <- FindMarkers(object = expr_test, ident.1 = 888, ident.2 = 999, test.use = "roc", only.pos = TRUE, min.pct = 0.25)
  }
  return(regionMk)
}


do_queryGene <- function(x) {
  x_inAnno <- x[x %in% geneID$symbol]
  x_inAnnoRmdup <- x[x %in% unique(geneID$symbol)]
  x_inExpr <- x[x %in% rownames(expr_data)]
  cat("Total:", length(x), "|", "In annotation:", length(x_inAnno), "|", "In rmdup annotation:", length(x_inAnnoRmdup), "|", "In expression table:", length(x_inExpr), "\n")
  cat("Unrecognized genes:", setdiff(x, x_inAnno), "\n")
}


DoHeatmap_new <- function (object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, 
          genes.use = NULL, genes.group = NULL, disp.min = -2.5, disp.max = 2.5, group.by = "ident", 
          switch_y = F, 
          group.order = NULL, draw.line = TRUE, col.low = "#FF00FF", 
          col.mid = "#000000", col.high = "#FFFF00", slim.col.label = FALSE, 
          remove.key = FALSE, rotate.key = FALSE, title = NULL, cex.col = 10, 
          cex.row = 10, group.label.loc = "bottom", group.label.rot = FALSE, 
          group.cex = 15, strip.text.y = 8, group.spacing = 0.15, panel.spacing.y = -0.2, assay.type = "RNA", 
          do.colBar = FALSE, colBar.y = 0.93, colBar.col = NULL, do.plot = TRUE, strip.text.x.top = 10, strip.text.y.display = F, strip.text.y.right = 4, legend.margin.for.colBar = margin(t = -12)) 
{
  if (is.null(x = data.use)) {
    if (use.scaled) {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "scale.data")
    }
    else {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "data")
    }
  }
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes_DF <- data.frame(row.names = genes.use, group = genes.group, stringsAsFactors = F)  # XXX
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  }
  else {
    cells.ident <- factor(x = FetchData(object = object, 
                                        cells.use = cells.use, vars.all = group.by)[, 1])
    names(x = cells.ident) <- cells.use
  }
  cells.ident <- factor(x = cells.ident, labels = intersect(x = levels(x = cells.ident), 
                                                            y = cells.ident))
  data.use <- data.use[genes.use, cells.use, drop = FALSE]
  if ((!use.scaled)) {
    data.use = as.matrix(x = data.use)
    if (disp.max == 2.5) 
      disp.max = 10
  }
  data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  data.use <- as.data.frame(x = t(x = data.use))
  data.use$cell <- rownames(x = data.use)
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  data.use <- data.use %>% melt(id.vars = "cell")
  names(x = data.use)[names(x = data.use) == "variable"] <- "gene"
  names(x = data.use)[names(x = data.use) == "value"] <- "expression"
  data.use$ident <- cells.ident[data.use$cell]
  if (!is.null(group.order)) {
    if (length(group.order) == length(levels(data.use$ident)) && 
        all(group.order %in% levels(data.use$ident))) {
      data.use$ident <- factor(data.use$ident, levels = group.order)
    }
    else {
      stop("Invalid group.order")
    }
  }
  breaks <- seq(from = min(data.use$expression), to = max(data.use$expression), 
                length = length(x = PurpleAndYellow()) + 1)
  data.use$gene <- with(data = data.use, expr = factor(x = gene, 
                                                       levels = rev(x = unique(x = data.use$gene))))
  data.use$cell <- with(data = data.use, expr = factor(x = cell, 
                                                       levels = cells.use))
  data.use$mkGroup <- genes_DF[as.character(data.use$gene), "group"]  # XXX
  data.use$mkGroup <- factor(data.use$mkGroup, levels = unique(genes.group))

  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  }
  else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  heatmap <- ggplot(data = data.use, mapping = aes(x = cell, 
                                                   y = gene, fill = expression)) + geom_tile() + scale_fill_gradient2(low = col.low, 
                                                                                                                      mid = col.mid, high = col.high, name = "Expression", 
                                                                                                                      guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) + 
    scale_y_discrete(position = "right" 
                     #labels = rev(genes.use)
                     ) + 
    theme(axis.line = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), strip.text.y = element_text(size = strip.text.y), 
          axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
          axis.title.x = element_blank())
  if (slim.col.label) {
    heatmap <- heatmap + theme(axis.title.x = element_blank(), 
                               axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                               axis.line = element_blank(), axis.title.y = element_blank(), 
                               axis.ticks.y = element_blank())
  }
  else {
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
  }
  if (!is.null(x = group.by)) {
    if (group.label.loc == "top") {
      switch <- NULL
    }
    else {
      switch <- "x"
    }
    if(switch_y) {  # XXX
      if(switch == "x") {
        switch <- "both"
      } else {
        switch <- "y"
      }
    }
    heatmap <- heatmap + facet_grid(facets = mkGroup~ident, drop = TRUE, # XXX
                                    space = "free", scales = "free", switch = switch 
    ) + scale_x_discrete(expand = c(0, 0), drop = TRUE)
    if (draw.line) {
      panel.spacing <- unit(x = group.spacing, units = "lines")
    }
    else {
      panel.spacing <- unit(x = 0, units = "lines")
    }
    heatmap <- heatmap + theme(strip.background = element_blank(), 
                               panel.spacing.x = panel.spacing, 
                               panel.spacing.y = unit(x = panel.spacing.y, units = "npc")
                               )
    if (group.label.rot) {
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
    }
  }
  if (remove.key) {
    heatmap <- heatmap + theme(legend.position = "none")
  }
  if (!is.null(x = title)) {
    heatmap <- heatmap + labs(title = title)
  }
  if (do.colBar) {
    ### pre
    if(group.label.rot) {
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = strip.text.x.top, r = 0, b = 0, l = 0)))
    } else {
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 0, margin = margin(t = strip.text.x.top, r = 0, b = 0, l = 0)))
    }
    heatmap <- heatmap + 
      theme(strip.text.y = element_text(margin = margin(r = strip.text.y.right))) + 
      theme(legend.position = "bottom", legend.justification = "center", legend.margin = legend.margin.for.colBar)
    if(! strip.text.y.display) {
      heatmap <- heatmap + theme(strip.text.y = element_blank())
    }
      #heatmap <- heatmap + theme()
    
    ###
    data.use_sub <- data.use[data.use$mkGroup == tail(levels(genes.group), 1), ]
    heatmap <- heatmap + geom_hline(data = data.use_sub, aes(yintercept = 1, color = ident), size = 2.5, show.legend = F)
    if(! is.null(colBar.col)) {
      heatmap <- heatmap + scale_color_manual(values = colBar.col)
    }
    ####
    g <- ggplotGrob(heatmap)
    
    # avoid clip
    for(i in which(grepl("strip-b", g$layout$name))){ g$grobs[[i]]$layout$clip <- "off" }
    
    groupNum <- length(unique(genes.group))
    botoom_ix <- grep("^panel-", g$layout$name)[1:length(grep("^panel-", g$layout$name)) %% groupNum == 0]
    
    for(i in seq_along(botoom_ix)) {
      layIx <- botoom_ix[i]
      objIx <- grep("GRID.segments", g$grobs[[layIx]]$childrenOrder)
      objID <- grep("GRID.segments", g$grobs[[layIx]]$childrenOrder, value = T)
      #cat(">", i, layIx, objIx, objID, "\n")
      
      # add new
      gt <- g$grobs[[grep(paste0("^strip-b-", i, "$"), g$layout$name)]]$grobs[[1]]
      gt$children <- gList(gt$children, added=g$grobs[[layIx]]$children[[objIx]])
      names(gt$children)[length(gt$children)] <- objID
      # change parameters
      gt$children[[length(gt$children)]]$y0 <- unit(colBar.y, units = "native")
      gt$children[[length(gt$children)]]$y1 <- unit(colBar.y, units = "native")
      gt$childrenOrder <- c(gt$childrenOrder, objID)
      g$grobs[[grep(paste0("^strip-b-", i, "$"), g$layout$name)]]$grobs[[1]] <- gt
      
      # remove old
      g$grobs[[layIx]]$children <- g$grobs[[layIx]]$children[-objIx]
      g$grobs[[layIx]]$childrenOrder <- g$grobs[[layIx]]$childrenOrder[-objIx]
    }
    
    heatmap <- ggplotify::as.ggplot(g)
    ####
  }
  if (do.plot) {
    heatmap
  }
  return(heatmap)
}
environment(DoHeatmap_new) <- asNamespace('Seurat')

DoHeatmap_object <- function (object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, 
                              genes.use = NULL, disp.min = -2.5, disp.max = 2.5, group.by = "ident", 
                              group.order = NULL, draw.line = TRUE, col.low = "#FF00FF", 
                              col.mid = "#000000", col.high = "#FFFF00", slim.col.label = FALSE, 
                              remove.key = FALSE, rotate.key = FALSE, title = NULL, cex.col = 10, 
                              cex.row = 10, group.label.loc = "bottom", group.label.rot = FALSE, 
                              group.cex = 15, group.spacing = 0.15, assay.type = "RNA", 
                              do.plot = TRUE) 
{
  if (is.null(x = data.use)) {
    if (use.scaled) {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "scale.data")
    }
    else {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "data")
    }
  }
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  }
  else {
    cells.ident <- factor(x = FetchData(object = object, 
                                        cells.use = cells.use, vars.all = group.by)[, 1])
    names(x = cells.ident) <- cells.use
  }
  cells.ident <- factor(x = cells.ident, labels = intersect(x = levels(x = cells.ident), 
                                                            y = cells.ident))
  data.use <- data.use[genes.use, cells.use, drop = FALSE]
  if ((!use.scaled)) {
    data.use = as.matrix(x = data.use)
    if (disp.max == 2.5) 
      disp.max = 10
  }
  data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  data.use <- as.data.frame(x = t(x = data.use))
  data.use$cell <- rownames(x = data.use)
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  data.use <- data.use %>% melt(id.vars = "cell")
  names(x = data.use)[names(x = data.use) == "variable"] <- "gene"
  names(x = data.use)[names(x = data.use) == "value"] <- "expression"
  data.use$ident <- cells.ident[data.use$cell]
  if (!is.null(group.order)) {
    if (length(group.order) == length(levels(data.use$ident)) && 
        all(group.order %in% levels(data.use$ident))) {
      data.use$ident <- factor(data.use$ident, levels = group.order)
    }
    else {
      stop("Invalid group.order")
    }
  }
  breaks <- seq(from = min(data.use$expression), to = max(data.use$expression), 
                length = length(x = PurpleAndYellow()) + 1)
  data.use$gene <- with(data = data.use, expr = factor(x = gene, 
                                                       levels = rev(x = unique(x = data.use$gene))))
  data.use$cell <- with(data = data.use, expr = factor(x = cell, 
                                                       levels = cells.use))
  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  }
  else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  
  ###
  return(data.use)
  ###
  
  # heatmap <- ggplot(data = data.use, mapping = aes(x = cell, 
  #                                                  y = gene, fill = expression)) + geom_tile() + scale_fill_gradient2(low = col.low, 
  #                                                                                                                     mid = col.mid, high = col.high, name = "Expression", 
  #                                                                                                                     guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) + 
  #   scale_y_discrete(position = "right", labels = rev(genes.use)) + 
  #   theme(axis.line = element_blank(), axis.title.y = element_blank(), 
  #         axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
  #         axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
  #         axis.title.x = element_blank())
  # if (slim.col.label) {
  #   heatmap <- heatmap + theme(axis.title.x = element_blank(), 
  #                              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
  #                              axis.line = element_blank(), axis.title.y = element_blank(), 
  #                              axis.ticks.y = element_blank())
  # }
  # else {
  #   heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
  # }
  # if (!is.null(x = group.by)) {
  #   if (group.label.loc == "top") {
  #     switch <- NULL
  #   }
  #   else {
  #     switch <- "x"
  #   }
  #   heatmap <- heatmap + facet_grid(facets = ~ident, drop = TRUE, 
  #                                   space = "free", scales = "free", switch = switch, 
  #   ) + scale_x_discrete(expand = c(0, 0), drop = TRUE)
  #   if (draw.line) {
  #     panel.spacing <- unit(x = group.spacing, units = "lines")
  #   }
  #   else {
  #     panel.spacing <- unit(x = 0, units = "lines")
  #   }
  #   heatmap <- heatmap + theme(strip.background = element_blank(), 
  #                              panel.spacing = panel.spacing)
  #   if (group.label.rot) {
  #     heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
  #   }
  # }
  # if (remove.key) {
  #   heatmap <- heatmap + theme(legend.position = "none")
  # }
  # if (!is.null(x = title)) {
  #   heatmap <- heatmap + labs(title = title)
  # }
  # if (do.plot) {
  #   heatmap
  # }
  # return(heatmap)
}
environment(DoHeatmap_object) <- asNamespace('Seurat')

DotPlot_new <- function (object, genes.plot, cols.use = c("lightgrey", "blue"), 
                         col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                         scale.by = "radius", scale.min = NA, scale.max = NA, group.by, 
                         plot.legend = FALSE, do.plot = TRUE, do.return = FALSE, x.lab.rot = FALSE, rev.x = F, rev.y = F, do.scale = T, breaks.use = NULL, limits.use = NULL, 
                         circle.color = "white") 
{
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  if(rev.y) {
    data.to.plot$id <- factor(data.to.plot$id, levels = rev(levels(data.to.plot$id)))
  }
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, 
                                                                            threshold = 0))
  if(do.scale) {
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
                                                                                 max = col.max, min = col.min))
  } else {
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp.scale = MinMax(data = avg.exp, 
                                  max = col.max, min = col.min))
  }
  
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                    levels = rev(x = genes.plot))
  if(rev.x) {
    data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, levels = rev(levels(data.to.plot$genes.plot)))
  }
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                 y = id)) + geom_point(mapping = aes(size = pct.exp, fill = avg.exp.scale), pch = 21, color = circle.color) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, 
                                                   scale.max)) + theme(axis.title.x = element_blank(), 
                                                                       axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_fill_distiller(palette = cols.use) # XXX
  }
  else {
    if(is.null(breaks.use)) {
      p <- p + scale_fill_gradient(low = cols.use[1], high = cols.use[2])  # XXX
    } else {
      p <- p + scale_fill_gradient(low = cols.use[1], high = cols.use[2], breaks = breaks.use, limits = limits.use)  # XXX
    }
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                              vjust = 0.5))
  }
  if (do.plot) {
    suppressWarnings(print(p))
  }
  if (do.return) {
    return(p)
  }
}
environment(DotPlot_new) <- asNamespace('Seurat')

SingleFeaturePlot_new <- function (data.use, feature, new.title = NULL, data.plot, pt.size, pch.use, cols.use, 
                                   dim.codes, min.cutoff, max.cutoff, coord.fixed, no.axes, 
                                   no.title = FALSE, no.legend, dark.theme, vector.friendly = FALSE, 
                                   png.file = NULL, png.arguments = c(10, 10, 100)) 
{
  if (vector.friendly) {
    previous_call <- blank_call <- png_call <- match.call()
    blank_call$pt.size <- -1
    blank_call$vector.friendly <- FALSE
    png_call$no.axes <- TRUE
    png_call$no.legend <- TRUE
    png_call$vector.friendly <- FALSE
    png_call$no.title <- TRUE
    blank_plot <- eval(blank_call, sys.frame(sys.parent()))
    png_plot <- eval(png_call, sys.frame(sys.parent()))
    png.file <- SetIfNull(x = png.file, default = paste0(tempfile(), 
                                                         ".png"))
    ggsave(filename = png.file, plot = png_plot, width = png.arguments[1], 
           height = png.arguments[2], dpi = png.arguments[3])
    to_return <- AugmentPlot(blank_plot, png.file)
    file.remove(png.file)
    return(to_return)
  }
  data.gene <- na.omit(object = data.frame(data.use[feature, 
                                                    ]))
  min.cutoff <- SetQuantile(cutoff = min.cutoff, data = data.gene)
  max.cutoff <- SetQuantile(cutoff = max.cutoff, data = data.gene)
  data.gene <- sapply(X = data.gene, FUN = function(x) {
    return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                  no = x))
  })
  data.gene <- sapply(X = data.gene, FUN = function(x) {
    return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                  no = x))
  })
  data.plot$gene <- data.gene
  if (length(x = cols.use) == 1) {
    brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
  }
  else {
    brewer.gran <- length(x = cols.use)
  }
  if (all(data.gene == 0)) {
    data.cut <- 0
  }
  else {
    data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.gene), 
                                                 breaks = brewer.gran)))
  }
  data.plot$col <- as.factor(x = data.cut)
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
  if (brewer.gran != 2) {
    if (length(x = cols.use) == 1) {
      p <- p + geom_point(mapping = aes(color = col), size = pt.size, 
                          shape = pch.use) + scale_color_brewer(palette = cols.use)
    }
    else {
      p <- p + geom_point(mapping = aes(color = col), size = pt.size, 
                          shape = pch.use) + scale_color_manual(values = cols.use)
    }
  }
  else {
    if (all(data.plot$gene == data.plot$gene[1])) {
      warning(paste0("All cells have the same value of ", 
                     feature, "."))
      p <- p + geom_point(color = cols.use[1], size = pt.size, 
                          shape = pch.use)
    }
    else {
      p <- p + geom_point(mapping = aes(color = gene), 
                          size = pt.size, shape = pch.use) + scale_color_gradientn(colors = cols.use, 
                                                                                   guide = guide_colorbar(title = feature))
    }
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  if (no.axes) {
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                   axis.text.y = element_blank(), axis.ticks = element_blank(), 
                   axis.title.x = element_blank(), axis.title.y = element_blank())
    if (!no.title) 
      ###
      if(is.null(new.title)) {
        p <- p + labs(title = feature, x = "", y = "")
      } else {
        p <- p + labs(title = new.title, x = "", y = "")
      }
      ###
    if (no.title) 
      p <- p + labs(x = "", y = "")
  }
  else {
    if (no.title) 
      p <- p + labs(x = dim.codes[1], y = dim.codes[2])
    if (!(no.title)) 
      ###
      if(is.null(new.title)) {
        p <- p + labs(title = feature, x = dim.codes[1], 
                      y = dim.codes[2])
      } else {
        p <- p + labs(title = new.title, x = dim.codes[1], 
                      y = dim.codes[2])
      }
      ###
  }
  if (no.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (coord.fixed) {
    p <- p + coord_fixed()
  }
  return(p)
}
environment(SingleFeaturePlot_new) <- asNamespace('Seurat')

FeaturePlot_new <- function (object, features.plot, new.title = NULL, min.cutoff = NA, max.cutoff = NA, 
          dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1, cols.use = c("yellow", 
                                                                            "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE, do.plot = T, 
          data.hover = "ident", do.identify = FALSE, reduction.use = "tsne", 
          use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE, 
          coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE, 
          vector.friendly = FALSE, png.file = NULL, png.arguments = c(10, 
                                                                      10, 100)) 
{
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = nCol)) {
    nCol <- 2
    if (length(x = features.plot) == 1) {
      nCol <- 1
    }
    if (length(x = features.plot) > 6) {
      nCol <- 3
    }
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
  }
  num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 
    1
  if (overlay | do.hover) {
    num.row <- 1
    nCol <- 1
  }
  par(mfrow = c(num.row, nCol))
  dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(GetCellEmbeddings(object = object, 
                                               reduction.type = reduction.use, dims.use = c(dim.1, 
                                                                                            dim.2), cells.use = cells.use))
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  data.plot$pt.size <- pt.size
  names(x = data.plot) <- c("x", "y")
  data.use <- t(x = FetchData(object = object, vars.all = features.plot, 
                              cells.use = cells.use, use.imputed = use.imputed))
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
                                                        ]), no = cutoff)
  }, cutoff = min.cutoff, feature = features.plot)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
                                                        ]), no = cutoff)
  }, cutoff = max.cutoff, feature = features.plot)
  check_lengths = unique(x = vapply(X = list(features.plot, 
                                             min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check_lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  if (overlay) {
    pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot, 
                            data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                            cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff, 
                            max.cutoff = max.cutoff, coord.fixed = coord.fixed, 
                            no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
  }
  else {
    if(is.null(new.title)) {
      new.title <- features.plot
    }
    pList <- mapply(FUN = SingleFeaturePlot_new, feature = features.plot, new.title = new.title, 
                    min.cutoff = min.cutoff, max.cutoff = max.cutoff, 
                    coord.fixed = coord.fixed, MoreArgs = list(data.use = data.use, 
                                                               data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                                                               cols.use = cols.use, dim.codes = dim.codes, 
                                                               no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme, 
                                                               vector.friendly = vector.friendly, png.file = png.file, 
                                                               png.arguments = png.arguments), SIMPLIFY = FALSE)
  }
  if (do.hover) {
    if (length(x = pList) != 1) {
      stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
    }
    if (is.null(x = data.hover)) {
      features.info <- NULL
    }
    else {
      features.info <- FetchData(object = object, vars.all = data.hover)
    }
    return(HoverLocator(plot = pList[[1]], data.plot = data.plot, 
                        features.info = features.info, dark.theme = dark.theme, 
                        title = features.plot))
  }
  else if (do.identify) {
    if (length(x = pList) != 1) {
      stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
    }
    return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, 
                          dark.theme = dark.theme))
  }
  else {
    if(do.plot) { # XXX
      print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
    }
  }
  ResetPar()
  if (do.return) {
    return(pList)
  }
}
environment(FeaturePlot_new) <- asNamespace('Seurat')
