
do_createDT <- function(expr_data, cellMetaData, do.norm = F, do.log = F, cell_num_cutoff = NULL) {
  if(! all(colnames(expr_data) == rownames(cellMetaData))) {
    stop("The cell sort in the two input are different, exit.")
  }
  if(do.norm) {
    cat("Normalize data...\n")
    expr_data <- sweep(expr_data, 2, colSums(expr_data), "/") * 1e6
  }
  if(do.log) {
    cat("Log transform data...\n")
    expr_data <- log2(expr_data + 1)
  }
  cellCluster_LS <- split(rownames(cellMetaData), cellMetaData$ident)
  # filtering
  if(! is.null(cell_num_cutoff)) {
    cellCluster_LS <- cellCluster_LS[lengths(cellCluster_LS) >= cell_num_cutoff]
  }
  
  res <- list()
  for(cell_type in names(cellCluster_LS)) {
    #cat(">", cell_type, "\n")
    res[[cell_type]][["data"]] <- expr_data[, colnames(expr_data) %in% cellCluster_LS[[cell_type]], drop = F]
  }
  return(res)
}

do_cor <- function(res_in, subsp = NULL, seed = 1, expr_cutoff = 0.1, mask = F, rm_outlier = F, method = "pearson", use = "all.obs") {
  cl <- makeCluster(min(length(res_in), 10), type = "FORK")
  expr_cor_LS <- parLapply(cl, res_in, function(x) {
    expr_input <- x[["data"]]
    # sub sampling
    if(! is.null(subsp)) {
      set.seed(seed)
      expr_input <- expr_input[, sample(1:ncol(expr_input), subsp)]
    }
    # filtering genes
    expr_input <- expr_input[rowSums(expr_input > 0) >= ncol(expr_input) * expr_cutoff, ]
    # mask zero
    if(mask) {
      expr_input[expr_input == 0] <- NA
    }
    # clac cor
    if(rm_outlier) {
      ngene <- nrow(expr_input)
      ncell <- ncol(expr_input)
      expr_cor <- matrix(NA, nrow = ngene, ncol = ngene)
      rownames(expr_cor) <- rownames(expr_input)
      colnames(expr_cor) <- rownames(expr_input)
      gene_mean <- rowMeans(expr_input)
      
      for(i in 1:ngene) {
        for(j in 1:ngene) {
          if(i==j) {
            expr_cor[i, j] <- 1
          } else {
            # identify outlier
            npro <- (expr_input[i, ] - gene_mean[i]) * (expr_input[j, ] - gene_mean[j])
            outlier <- which.max(npro)
            # recalculate mean without outlier
            if(method == "bicor") {
              expr_cor[i, j] <- WGCNA::cor(expr_input[i, -outlier], expr_input[j, -outlier], use = use)
            } else {
              expr_cor[i, j] <- cor(expr_input[i, -outlier], expr_input[j, -outlier], method = method, use = use)
            }
          }
        }
      }
    } else {
      if(method == "bicor") {
        expr_cor <- WGCNA::bicor(t(expr_input), use = use)
      } else {
        expr_cor <- cor(t(expr_input), method = method, use = use)
      }
    }
    return(expr_cor)
  })
  stopCluster(cl); rm(cl)
  for(cell_type in names(res_in)) {
    res_in[[cell_type]][["cor"]] <- expr_cor_LS[[cell_type]]
  }
  return(res_in)
}

do_hc <- function(res_in, use.abs = F, method = "average", rm.cor = T) {
  cl <- makeCluster(min(length(res_in), 10), type = "FORK")
  hc_LS <- parLapply(cl, names(res_in), function(x) {
    print(x)
    expr_cor <- res_in[[x]][["cor"]]
    if(use.abs) {
      ds <- as.dist(1 - abs(expr_cor))
    } else {
      ds <- as.dist(1 - expr_cor)
    }
    hc <- hclust(ds, method = method)
    expr_cor_clustered <- expr_cor[hc$order, hc$order]
    y <- list(hc = hc, cor_cld = expr_cor_clustered)
    return(y)
  })
  stopCluster(cl); rm(cl)
  names(hc_LS) <- names(res_in)
  for(cell_type in names(hc_LS)) {
    res_in[[cell_type]][["hc"]] <- hc_LS[[cell_type]][["hc"]]
    res_in[[cell_type]][["cor_cld"]] <- hc_LS[[cell_type]][["cor_cld"]]
    if(rm.cor) {
      res_in[[cell_type]][["cor"]] <- NULL
    }
  }
  return(res_in)
}

do_diffCut <- function(hc_in, con_cutoff = 2.5, len_cutoff = 1) {
  hc_DF <- data.frame(ka = hc_in$merge[, 1], kb = hc_in$merge[, 2], height = hc_in$height)
  hc_DF$top_len <- sapply(1:nrow(hc_DF), function(i) {
    ix <- c(which(hc_DF[, 1] == i), which(hc_DF[, 2] == i))
    if(length(ix) > 0) {
      y <- hc_DF[ix, 3] - hc_DF[i, 3]
    } else {
      y <- NA
    }
    return(y)
  })
  
  hc_DF$top_cut <- hc_DF$height <= con_cutoff & hc_DF$top_len >= len_cutoff & hc_DF$ka > 0 & hc_DF$kb > 0
  mgid_cut <- which(hc_DF$top_cut)
  
  # get elements in each mg
  hc_LS <- list()
  for(i in 1:nrow(hc_DF)) {
    #cat(">", i, "\n")
    xs <- as.numeric(hc_DF[i, 1:2])
    xe <- xs[xs < 0]  # element
    xc <- xs[xs > 0]  # cluster
    for(xi in xc) {
      #cat(">>", xi, "\n")
      xe <- c(xe, hc_LS[[xi]])
    }
    hc_LS[[i]] <- xe
  }
  
  cut_out <- cutree(hc_in, h = con_cutoff)
  for(i in mgid_cut) {
    cut_out[ hc_in$labels[-hc_LS[[i]]] ] <- paste(cut_out[ hc_in$labels[-hc_LS[[i]]] ], i, sep = ".")
  }
  
  return(cut_out)
}

do_detectModule <- function(res_in, method = "constant", h_cutoff, size_cutoff = 10, avgCor_cutoff = NULL, write.output = T) {
  # read gene type info
  gene_type <- read.table("/rd/user/tianf/06-Human_cell_atlas/Genomes/human/gene_type_class.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 2)
  colnames(gene_type) <- c("ensembl_id", "type", "class")
  gene_type$class <- factor(gene_type$class, levels = unique(gene_type$class)[c(4,3,2,1,5)])
  # read tf list
  tf_list <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/AnimalTFDB/human_TF_list.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]
  
  for(cell_type in names(res_in)) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    
    hc <- res_in[[cell_type]][["hc"]]
    expr_cor_clustered <- res_in[[cell_type]][["cor_cld"]]
    # cut tree
    if(method == "constant") {
      cut_res <- cutree(tree = hc, h = h_cutoff)
    } else if (method == "dynamic") {
      cut_res <- cutreeDynamic(dendro = hc, cutHeight = h_cutoff, minClusterSize = size_cutoff, method = "tree", deepSplit = F)
    } else {
      stop("The method should be either \"constant\" or \"dynamic\".")
    }
    cluster_table <- data.frame(gene = hc$labels, cluster = cut_res, stringsAsFactors = F)
    cluster_table$cluster <- paste0("M", cluster_table$cluster)
    cluster_table[cluster_table$cluster == "M0", "cluster"] <- paste0("M0_", cluster_table[cluster_table$cluster == "M0", "gene"])
    cluster_gene <- split(cluster_table$gene, cluster_table$cluster)[unique(cluster_table$cluster)]
    # add size
    cluster_size <- data.frame(size = table(cluster_table$cluster))
    colnames(cluster_size) <- c("cluster", "size")
    cluster_table <- merge(cluster_table, cluster_size, by.x = 2, by.y = 1, sort = F)
    # add boundary
    cluster_table <- merge(data.frame(gene=rownames(expr_cor_clustered), stringsAsFactors = F), cluster_table, by.x = 1, by.y = 2, sort = F) # sort by gene ID in clustered cor
    cluster_bdl <- c(1, cumsum(cluster_table[! duplicated(cluster_table$cluster), "size"])[-length(unique(cluster_table$cluster))] + 1)
    cluster_bdr <- cumsum(cluster_table[! duplicated(cluster_table$cluster), "size"])
    cluster_bd <- data.frame(cluster=unique(cluster_table$cluster), boundary=paste(cluster_bdl, cluster_bdr, sep = ":"), stringsAsFactors = F)
    cluster_table <- merge(cluster_table, cluster_bd, by = "cluster", sort = F)
    # add expr stat
    cluster_exprStat <- t(sapply(cluster_gene, function(x) { y0 <- res_in[[cell_type]]$data[x, ]; y <- c(mean(unlist(y0)), mean(apply(y0, 1, sd) / rowMeans(y0))); return(y) }))
    colnames(cluster_exprStat) <- c("avgExpr", "cvExpr")
    cluster_table <- merge(cluster_table, cluster_exprStat, by.x = "cluster", by.y = 0, sort = F)
    # add cor stat
    cluster_corStat <- t(sapply(cluster_gene, function(x) { y0 <- expr_cor_clustered[x, x]; y1 <- y0[upper.tri(y0)]; y <- c(mean(y1), sd(y1) / mean(y1)); return(y) }))
    colnames(cluster_corStat) <- c("avgCor", "cvCor")
    cluster_table <- merge(cluster_table, cluster_corStat, by.x = "cluster", by.y = 0, sort = F)
    # add gene type stat
    cluster_gtStat <- t(sapply(cluster_gene, function(x) { y <- table(gene_type[x, "class"]) / length(x); return(y) }))
    cluster_table <- merge(cluster_table, cluster_gtStat, by.x = "cluster", by.y = 0, sort = F)
    # add tf stat
    cluster_tfStat <- data.frame(tf = sapply(cluster_gene, function(x) { y <- sum(x %in% tf_list) / length(x); return(y) }))
    cluster_table <- merge(cluster_table, cluster_tfStat, by.x = "cluster", by.y = 0, sort = F)
    # filtering
    cluster_table[grep("M0_", cluster_table$cluster), "cluster"] <- "M0"
    cat("Raw:", length(unique(cluster_table$cluster)), "modules with", nrow(cluster_table), "genes\n")
    cluster_table_ftd <- subset(cluster_table, cluster != "M0" & size >= size_cutoff)
    if(! is.null(avgCor_cutoff)) {
      cluster_table_ftd <- subset(cluster_table_ftd, avgCor >= avgCor_cutoff)
    }
    cat("Clean:", length(unique(cluster_table_ftd$cluster)), "modules with", nrow(cluster_table_ftd), "genes\n")
    if(write.output) {
      write.table(x = cluster_table_ftd, file = paste0("03-expression/merged/geneModule/geneModule_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    }
    res_in[[cell_type]][["cl_table_ftd"]] <- cluster_table_ftd
    # filtering cor
    #expr_cor_ftd <- expr_cor_clustered[rownames(expr_cor_clustered)%in%cluster_table_ftd$gene, colnames(expr_cor_clustered)%in%cluster_table_ftd$gene]
    #write.table(x = expr_cor_ftd, file = paste0(OUT, "/geneCorrelation_", cell_type_label, ".txt", row.names = T, col.names = NA, quote = F, sep = "\t")
    #res_in[[cell_type]][["cor_cld_ftd"]] <- expr_cor_ftd
    # module list
    cluster_list_ftd <- unique(cluster_table_ftd[, -2])
    rownames(cluster_list_ftd) <- NULL
    res_in[[cell_type]][["cl_list_ftd"]] <- cluster_list_ftd
  }
  return(res_in)
}

do_plotCoordCor <- function(res_in, ctype1, ctype2, do_plot = F, sub_num = NULL, show_rownames = F, show_colnames = F) {
  gene1 <- res_in[[ctype1]]$cl_table_ftd[, "gene"]
  gene2 <- res_in[[ctype2]]$cl_table_ftd[, "gene"]
  genes <- intersect(gene1, gene2)
  cat("In-module gene number:", "\n")
  cat("  Total:", length(gene1), "from", ctype1, ",", length(gene2), "from", ctype2, "\n")
  cat("  Overlapped:", length(genes), "\n")
  cat("  Exclusive:", length(setdiff(gene1, gene2)), length(setdiff(gene2, gene1)), "\n")
  
  # overlapping stat using barplot
  if(do_plot) {
    intersect_DF <- data.frame(ctype=c(ctype1, "overlapped", ctype2), number=c(length(gene1)-length(genes), length(genes), length(gene2)-length(genes)))
    intersect_DF$ctype <- factor(intersect_DF$ctype, levels = c(ctype2, "overlapped", ctype1))  # should rev
    rect_DF <- data.frame(x = c(0.4, 0.4, 1.6, 1.6), 
                          number = c(intersect_DF$number[1], sum(intersect_DF$number[1:2]), sum(intersect_DF$number[1:2]), intersect_DF$number[1]), 
                          ctype = NA, group = 1:4)
    
    p0 <- ggplot(intersect_DF, aes(x = 1, y = number, fill = ctype)) + geom_bar(stat='identity', show.legend = F) + coord_flip() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + 
      scale_x_continuous(expand = c(0.5, 0)) + scale_y_continuous(expand = c(0.1, 0)) + 
      ylab("Gene number") + 
      geom_polygon(data = rect_DF, aes(x = x), color = "black", linetype = "dashed", fill = NA)
    
    cor_1 <- res_in[[ctype1]]$cor_cld_ftd[genes, genes]
    cor_2 <- res_in[[ctype2]]$cor_cld_ftd[genes, genes]
    if(! is.null(sub_num)) {
      warning("Only a fraction of genes were used for plot.")
      cor_1 <- cor_1[sub_num, sub_num]
      cor_2 <- cor_2[sub_num, sub_num]
    }
    p1 <- pheatmap(mat = cor_1, cluster_rows = F, cluster_cols = F, breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1),
                   color = c("blue",colorRampPalette(c("blue", "white", "red"))(100),"red"),
                   show_rownames = show_rownames, show_colnames = show_colnames, legend = F,
                   main = ctype1, silent = T)
    p2 <- pheatmap(mat = cor_2, cluster_rows = F, cluster_cols = F, breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1),
                   color = c("blue",colorRampPalette(c("blue", "white", "red"))(100),"red"),
                   show_rownames = show_rownames, show_colnames = show_colnames, legend = F,
                   main = ctype2, silent = T)
    #png(paste0("03-expression/merged/geneModule/", samplingPos, "/coordCor_", ctype1, "_", ctype2, ".png"), height = 1200, width = 1600, res = 300)
    #pdf(paste0("03-expression/merged/geneModule/", samplingPos, "/coordCor_", ctype1, "_", ctype2, ".pdf"), height = 5.5, width = 8)
    #print(p0)
    gridExtra::grid.arrange(p0, p1[[4]], p2[[4]], layout_matrix = rbind(c(1,1), c(2,3), c(2,3), c(2,3), c(2,3)))
    #dev.off()
  }
}

expandCorMat <- function(cand, ref, value = NA, fill_diag = T) {
  # expand cand according ref
  # remove inconsistent names in cand
  cand <- cand[rownames(cand) %in% rownames(ref), colnames(cand) %in% colnames(ref)]
  # create the full matrix
  ko <- ref
  ko[,] <- value
  # fill the full matrix
  ind_row <- match(rownames(cand), rownames(ko))
  ind_col <- match(colnames(cand), colnames(ko))
  ko[ind_row, ind_col] <- cand
  if(fill_diag) {
    diag(ko) <- 1
  }
  return(ko)
}

corMatStat <- function(mt_in, cond, FUN = function(x) { mean(x, na.rm = T)}) {
  # stat for each x~y defined region
  mt <- mt_in
  diag(mt) <- NA
  # rows
  corm <- aggregate(mt, by = cond, FUN)
  rownames(corm) <- corm[, 1]
  corm <- corm[, -1]
  # cols
  corm <- aggregate(x = t(corm), by = cond, FUN)
  rownames(corm) <- corm[, 1]
  corm <- corm[, -1]
  # final
  corm <- t(corm)
  return(corm)
}

diffMd <- function(cor1, cor2, ct1, cutoff = 0.05) {
  ct1$cluster <- factor(ct1$cluster, levels = unique(ct1$cluster))
  ct1_LS <- split(ct1$gene, ct1$cluster)
  out <- sapply(ct1_LS, function(x) {
    cob1 <- cor1[x, x]
    cob2 <- cor2[x, x]
    coa1 <- cob1[upper.tri(cob1)]
    coa2 <- cob2[upper.tri(cob2)]
    avg1 <- mean(coa1, na.rm = T)
    avg2 <- mean(coa2, na.rm = T)
    fc <- avg2 / avg1
    pvalue <- wilcox.test(coa1, coa2, paired = T)$p.value
    return(c(avg1, avg2, fc, pvalue))
  })
  out <- data.frame(t(out))
  rownames(out) <- names(ct1_LS)
  colnames(out) <- c("avg_1", "avg_2", "fc", "pvalue")
  out$qvalue <- p.adjust(out$pvalue, method = "BH")
  out$holdOn <- ifelse(! (out$fc < 1 & out$qvalue < 0.05), 1, 0)
  return(out)
}

do_cmpMd <- function(res_in, ctype1, ctype2, do_plot = F, fontsize_row = 6, vmin = -0.3, vmax = 0.3) {
  cor1 <- res_in[[ctype1]][["cor_cld"]]
  cor2 <- res_in[[ctype2]][["cor_cld"]]
  cor2 <- expandCorMat(cor2, cor1)  # expand
  ct1 <- res_in[[ctype1]][["cl_table_ftd"]]
  # diff module
  diff_DF <- diffMd(cor1, cor2, ct1, cutoff = 0.1)
  colnames(diff_DF)[1:2] <- c(ctype1, ctype2)
  cat(sum(diff_DF$holdOn), "/", length(diff_DF$holdOn), "modules are still correlated.\n")
  
  if(do_plot) {
    #pdf(paste0("03-expression/merged/geneModule/", samplingPos, "/diffMd_", ctype1, "_", ctype2, ".pdf"), height = 5.5, width = 3)
    anno_row_DF <- data.frame(diff_DF[, "holdOn", drop = F])
    anno_row_DF$holdOn <- factor(anno_row_DF$holdOn)
    pheatmap_new(diff_DF[, 1:2], cluster_rows = F, cluster_cols = F, 
                 breaks = c(-1, seq(from = vmin, to = vmax, length.out = 100), 1), 
                 color = c("blue",colorRampPalette(c("blue", "white", "red"))(100), "red"), 
                 annotation_row = anno_row_DF, annotation_colors = list(holdOn = c(`0` = "grey", `1` = "green")), 
                 annotation_names_row = F, angle_col = 45, fontsize_row = fontsize_row)
    #dev.off()
  }
  return(diff_DF)
}

do_plotCorHeatmap <- function(res_in, ctype1, ctype2 = NULL, mdid = NULL, mp_in = NULL, mpid = NULL, mgid = NULL, 
                              do_plot = T, rm.upper = F, na_col = "#DDDDDD", vmin = -0.6, vmax = 0.6, do_highlights = F, 
                              main = NA, show_rownames = F, show_colnames = F, fontsize_row = 6, fontsize_col = 6, show.legend = T, 
                              do.print = T, do.return = T) {
  cor1 <- res_in[[ctype1]][["cor_cld"]]
  if(! is.null(mdid)) {
    ct1 <- res_in[[ctype1]][["cl_table_ftd"]]
    ct1_sub <- subset(ct1, cluster %in% mdid, drop = F)
    ct1_sub <- do.call("rbind", split(ct1_sub, ct1_sub$cluster)[mdid])  # sort by mdid
    genes <- ct1_sub$gene
  } else if(! is.null(mpid)) {
    genes <- mp_in[[mpid]]
  } else if(! is.null(mgid)) {
    genes <- mgid
  } else {
    stop("Either mdid or mpList + mpid should be specified, exit.")
  }
  if(is.null(genes)) {
    stop("The mdid or mpid is missing, stop.")
  }
  cor_comb <- matrix(NA, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
  cor_sub1 <- expandCorMat(cand = cor1, ref = cor_comb) # expand
  cor_comb <- cor_sub1
  if(! is.null(ctype2)) {
    cor2 <- res_in[[ctype2]][["cor_cld"]]
    cor_sub2 <- expandCorMat(cand = cor2, ref = cor_comb)  # expand
    cor_comb[upper.tri(cor_comb)] <- cor_sub2[upper.tri(cor_sub2)]
  }
  if(do_plot) {
    if(do_highlights) {
      ct1_sub_uniq <- unique(ct1_sub[, c("cluster", "size")])
      ind1 <- c(1, cumsum(ct1_sub_uniq$size)[-nrow(ct1_sub_uniq)] + 1)
      ind2 <- cumsum(ct1_sub_uniq$size)
      ind_highlights <- data.frame(st = ind1, ed = ind2, stringsAsFactors = F)
    } else {
      ind_highlights <- NULL
    }
    # limit value space
    cor_comb[cor_comb < vmin] <- vmin
    cor_comb[cor_comb > vmax] <- vmax
    # tri plot
    if(rm.upper) {
      cor_comb[upper.tri(cor_comb)] <- 0
    }
    
    cor_comb_melted <- melt(cor_comb)
    cor_comb_melted$Var1 <- factor(cor_comb_melted$Var1, levels = rev(unique(cor_comb_melted$Var1)))
    p <- ggplot(cor_comb_melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile(show.legend = show.legend) + 
      scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100), limits = c(vmin, vmax), breaks = c(vmin, vmin / 2, 0, vmax / 2, vmax), na.value = na_col) + 
      scale_x_discrete(position = "top") + labs(fill = "Correlation") + xlab(ctype2) + ylab(ctype1) + ggtitle(main) + 
      theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
      theme(legend.margin = margin(l = -5), aspect.ratio = 1)
    if(! is.null(ind_highlights)) {
      p <- p + geom_rect(data = ind_highlights, aes(xmin = st - 0.5, xmax = ed + 0.5, ymin = length(genes) + 1 - (st - 0.5), ymax = length(genes) + 1 - (ed + 0.5)), fill = NA, color = "black", inherit.aes = F)
    }
    if(is.null(ctype2)) {
      p <- p + ylab(NULL)
    }
    if(show_rownames) {
      p <- p + theme(axis.text.y = element_text(size = fontsize_row))
    }
    if(show_colnames) {
      p <- p + theme(axis.text.x = element_text(size = fontsize_col))
    }
    if(is.na(main)) {
      p <- p + ggtitle(NULL)
    }
    if(do.print) {
      print(p)
    } else if(do.return) {
      return(p)
    }
  }
}

do_plotExprHeatmap <- function(res_in, ctype1, ctype2 = NULL, mdid = NULL, mp_in, mpid = NULL, subsp = 500, seed = 1, mask0 = F, do.log = T, scale = "row", 
                               vmin = -2.5, vmax = 2.5, 
                               main = NA, show_rownames = F, show_colnames = F, fontsize_row = 6, fontsize_col = 6, show.legend = T, 
                               do.print = T, do.return = T, short_name = F, asp.ratio = 2.5) {
  if(! is.null(mdid)) {
    ct1 <- res_in[[ctype1]][["cl_table_ftd"]]
    ct1_sub <- subset(ct1, cluster %in% mdid, drop = F)
    ct1_sub <- do.call("rbind", split(ct1_sub, ct1_sub$cluster)[mdid])  # sort by mdid
    genes <- ct1_sub$gene
  } else if(! is.null(mpid)) {
    genes <- mp_in[[mpid]]
  } else {
    stop("Either mdid or mpid should be specified, exit.")
  }
  if(is.null(genes)) {
    stop("The mdid or mpid is missing, stop.")
  }
  expr1 <- res_in[[ctype1]][["data"]]
  set.seed(seed)
  expr_sub1 <- expr1[genes, sample(1:ncol(expr1), subsp)]
  expr_sub1 <- expr_sub1[, order(colMeans(sweep(expr_sub1, 1, rowSums(expr_sub1), "/")))]
  expr_comb <- expr_sub1
  if(! is.null(ctype2)) {
    expr2 <- res_in[[ctype2]][["data"]]
    set.seed(seed)
    expr_sub2 <- expr2[genes, sample(1:ncol(expr2), subsp)]
    expr_sub2 <- expr_sub2[, order(colMeans(sweep(expr_sub2, 1, rowSums(expr_sub2), "/")))]
    expr_comb <- cbind(expr_sub1, expr_sub2)
  }
  if(mask0) {
    expr_comb[expr_comb == 0] <- NA
  }
  if(do.log) {
    expr_comb <- log2(expr_comb + 1)
  }
  # scale
  expr_comb <- switch(scale, none = expr_comb, row = pheatmap:::scale_rows(expr_comb), column = t(pheatmap:::scale_rows(t(expr_comb))))
  # limit value space
  expr_comb[expr_comb < vmin] <- vmin
  expr_comb[expr_comb > vmax] <- vmax
  expr_comb_melted <- melt(as.matrix(expr_comb))
  expr_comb_melted$Var1 <- factor(expr_comb_melted$Var1, levels = rev(unique(expr_comb_melted$Var1)))
  if(! is.null(ctype2)) {
    expr_comb_melted$ct <- c(rep(ctype1, nrow(expr_sub1) * ncol(expr_sub1)), rep(ctype2, nrow(expr_sub2) * ncol(expr_sub2)))
    expr_comb_melted$ct <- factor(expr_comb_melted$ct, levels = c(ctype1, ctype2))
    if(short_name) {
      split_str <- strsplit(levels(expr_comb_melted$ct), split = ".", fixed = T)
      split_str_ts <- sapply(split_str, "[", 1)
      split_str_ct <- sapply(split_str, "[", 2)
      levels(expr_comb_melted$ct) <- paste(split_str_ts, substr(split_str_ct, 1, 1), sep = ".")
    }
  }
  p <- ggplot(expr_comb_melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile(show.legend = show.legend) + 
    scale_fill_gradientn(colours = colorRampPalette(c("purple3", "white","palegreen3"))(100), na.value = "grey80", limits = c(vmin, vmax), breaks = c(vmin, 0, vmax)) + 
    scale_x_discrete(position = "top") + labs(fill = "Expression") + xlab(NULL) + ylab(NULL) + ggtitle(main) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
    theme(plot.margin = margin(6.5,7,4,7)) + 
    theme(aspect.ratio = asp.ratio)
  if(! is.null(ctype2)) {
    p <- p + facet_grid(. ~ ct, scales = "free_x", switch = "x") + theme(strip.background = element_blank(), strip.text.x = element_text(margin = margin(4.4, 4.4, 4.4, 4.4)))
  }
  if(show_rownames) {
    p <- p + theme(axis.text.y = element_text(size = fontsize_row))
  }
  if(show_colnames) {
    p <- p + theme(axis.text.x = element_text(size = fontsize_col))
  }
  if(is.na(main)) {
    p <- p + ggtitle(NULL)
  }
  if(do.print) {
    print(p)
  } else if(do.return) {
    return(p)
  }
}

do_plotHclust <- function(res_in, ctype1, ctype2 = NULL, mdid = NULL, method = "average", 
                          horiz = F, do_rev = F, rm_labels = F, main = NULL, xlab = NA, sub = NA, do_return = F) {
  if(is.null(mdid)) {
    hc_sub1 <- res_in[[ctype1]]$hc
  } else {
    genes <- subset(res_in[[ctype1]]$cl_table_ftd, cluster %in% mdid, "gene", drop = T)
    cor_sub1 <- res_in[[ctype1]]$cor_cld[genes, genes]
    hc_sub1 <- hclust(as.dist(1-cor_sub1), method = method)
  }
  if(do_rev) {
    dd_sub1 <- as.dendrogram(hc_sub1)
    hc_sub1 <- as.hclust(rev(dd_sub1))
  }
  if(horiz) {
    par(mar=c(3,1,4,5))
  } else if(rm_labels) {
    par(mar=c(1,4,4,2))
  } else {
    par(mar=c(8,4,4,2))
  }
  if(rm_labels) {
    plot(hc_sub1, main = main, labels = F, xlab = xlab, sub = sub)
  } else {
    plot(hc_sub1, main = main, xlab = xlab, sub = sub)
  }
  par(mar=c(5,4,4,2) + 0.1)
  if(do_return) { return(hc_sub1) }
}

diffGene <- function(expr_data, cellMetaData, ctype1, ctype2) {
  if(any(colnames(expr_data) != rownames(cellMetaData))) {
    exit("Inconsistent cell names between expression data and meta table.")
  }
  out <- apply(expr_data, 1, function(x) {
    x1 <- x[cellMetaData$expr.ident==ctype1]
    x2 <- x[cellMetaData$expr.ident==ctype2]
    pvalue <- wilcox.test(x1, x2)$p.value
    logFC <- log10(mean(x2) / mean(x1))
    return(c(pvalue, logFC))
  })
  out <- data.frame(t(out))
  colnames(out) <- c("pvalue", "logFC")
  out$qvalue <- p.adjust(out$pvalue, method = "BH")
  return(out)
}

geneID_mapping <- function() {
  # get ID table from biomart
  library("biomaRt")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  id_DF <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), mart=mart)
  #id_DF_rmdup <- id_DF[! duplicated(id_DF$ensembl_gene_id), ]
  
  # read gencode gene list
  genes <- read.table("/rd/user/tianf/06-Human_cell_atlas/Genomes/human/gene_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)
  colnames(genes) <- c("gencode", "symbol")
  genes$ensemblGene <- sapply(strsplit(genes$gencode, ".", fixed=T), "[", 1)
  #genes <- genes[! genes$symbol %in% genes[duplicated(genes$symbol), "symbol"], ]
  # merge
  genes_DF <- merge(x = genes, y = id_DF, by.x = "ensemblGene", by.y = "ensembl_gene_id", sort = F, all.x = T)
  #genes_DF <- genes_DF[match(genes$gencode, genes_DF$gencode), ]
  #all(genes$gencode == genes_DF$gencode)
  write.table(x = genes_DF, file = "/rd/user/tianf/06-Human_cell_atlas/Data/human_id_table.txt", row.names = F, col.names = T, quote = F, sep = "\t")
}

enricher_new <- function (gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, 
                          minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE, 
                          TERM2NAME = NA, user_data = NULL) 
{
  if(is.null(user_data)) {
    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
  } else {
    USER_DATA <- user_data
  }
  enricher_internal(gene = gene, pvalueCutoff = pvalueCutoff, 
                    pAdjustMethod = pAdjustMethod, universe = universe, 
                    minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, 
                    USER_DATA = USER_DATA)
}
environment(enricher_new) <- asNamespace('clusterProfiler')

do_prepGO <- function() {
  go_anno <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/GO/goa_human.gaf", header = F, sep = "\t", stringsAsFactors = F, 
                        comment.char = "!", quote = "")
  dim(go_anno)
  # rmdup
  go_anno <- go_anno[, c(3,5,9)]
  go_anno <- go_anno[! duplicated(go_anno), ]
  # get ancestor node
  go_ancestor <- GOfuncR::get_parent_nodes(go_ids = unique(go_anno$V5))
  # check missed GO ID
  length(unique(go_anno$V5))
  length(unique(go_ancestor$child_go_id))
  setdiff(unique(go_anno$V5), unique(go_ancestor$child_go_id))
  # anno with ancestor
  go_anno_an <- merge(go_anno, go_ancestor, by.x = "V5", by.y = "child_go_id", sort = F)
  go_anno_an <- go_anno_an[, 2:5]
  go_anno_an <- go_anno_an[! duplicated(go_anno_an), ]
  GO2gene <- go_anno_an[, c(3,1,2)]
  colnames(GO2gene) <- c("GO_ID", "Gene_name", "aspect")
  GO2name <- go_anno_an[, c(3,4,2)]
  GO2name <- GO2name[! duplicated(GO2name), ]
  colnames(GO2name) <- c("GO_ID", "GO_name", "aspect")
  # build anno
  for(i in c("P", "F", "C")) {
    cat(">", i, "\n")
    GO2gene_sub <- subset(GO2gene, aspect == i)[, 1:2]
    GO2name_sub <- subset(GO2name, aspect == i)[, 1:2]
    anno_build <- clusterProfiler:::build_Anno(path2gene = GO2gene_sub, path2name = GO2name_sub)
    saveRDS(object = anno_build, file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/GO/GO_", i, ".rds"))
  }
}

do_enrichGO <- function(res_in, ctype = NULL, mdid = NULL, aspect = "P", ncpu = 1) {
  # read GO anno
  gs <- readRDS(file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/GO/GO_", aspect, ".rds"))
  
  if(! is.null(ctype)) {
    cell_type_list <- ctype
  } else {
    cell_type_list <- names(res_in)
  }
  
  cl <- makeCluster(ncpu, type = "FORK")
  for(cell_type in cell_type_list) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    cluster_table_ftd <- res_in[[cell_type]]$cl_table_ftd
    cluster_table_ftd$cluster <- factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))
    cluster_gene <- split(cluster_table_ftd$gene, cluster_table_ftd$cluster)
    if(! is.null(mdid)) {
      cluster_gene <- cluster_gene[mdid]
    }
    enriched_LS <- parLapply(cl, seq_along(cluster_gene), function(i) {
      cat(">>", i, "\n")
      geneInput <- cluster_gene[[i]]
      y0 <- enricher_new(gene = geneInput, TERM2GENE = NA, TERM2NAME = NA, user_data = gs, 
                         pAdjustMethod = "BH", minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 0.05)
      y0 <- as.data.frame(y0)
      y0 <- data.frame(cluster = rep(names(cluster_gene)[i], nrow(y0)), y0, stringsAsFactors = F)
      # filtering
      y <- subset(y0, Count >= 3 & p.adjust <= 0.05)
      return(y)
    })
    enrich_res <- do.call("rbind", enriched_LS)
    write.table(x = enrich_res, file = paste0(OUT, "/enrich_GO_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    res_in[[cell_type]][["enrich_GO"]] <- enrich_res
    # add stat
    res_in[[cell_type]]$cl_table_ftd$enrich_GO <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res$cluster
    res_in[[cell_type]]$cl_list_ftd$enrich_GO <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res$cluster
  }
  stopCluster(cl); rm(cl)
  return(res_in)
}

do_prepKEGG <- function() {
  term2gene <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/KEGG/KEGG_pathway2gene_formatted.txt", header = F, sep = "\t", stringsAsFactors = F)
  dim(term2gene)
  colnames(term2gene) <- c("pathway", "gene")
  term2gene$pathway <- gsub("path:", "", term2gene$pathway)
  
  term2name <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/KEGG/KEGG_pathway2name.txt", header = F, sep = "\t", stringsAsFactors = F)
  dim(term2name)
  colnames(term2name) <- c("pathway", "name")
  term2name$pathway <- gsub("path:", "", term2name$pathway)
  term2name$name <- gsub(" - Homo sapiens (human)", "", term2name$name, fixed = T)
  
  # build anno
  anno_build <- clusterProfiler:::build_Anno(path2gene = term2gene, path2name = term2name)
  saveRDS(object = anno_build, file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/KEGG/KEGG_", "pathway", ".rds"))
}

do_enrichKEGG <- function(res_in, ctype = NULL, mdid = NULL, ncpu = 1) {
  # read anno
  gs <- readRDS(file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/KEGG/KEGG_", "pathway", ".rds"))
  
  if(! is.null(ctype)) {
    cell_type_list <- ctype
  } else {
    cell_type_list <- names(res_in)
  }
  
  cl <- makeCluster(ncpu, type = "FORK")
  for(cell_type in cell_type_list) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    cluster_table_ftd <- res_in[[cell_type]]$cl_table_ftd
    cluster_table_ftd$cluster <- factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))
    cluster_gene <- split(cluster_table_ftd$gene, cluster_table_ftd$cluster)
    if(! is.null(mdid)) {
      cluster_gene <- cluster_gene[mdid]
    }
    enriched_LS <- parLapply(cl, seq_along(cluster_gene), function(i) {
      cat(">>", i, "\n")
      geneInput <- cluster_gene[[i]]
      y0 <- enricher_new(gene = geneInput, TERM2GENE = NA, TERM2NAME = NA, user_data = gs, 
                         pAdjustMethod = "BH", minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 0.05)
      y0 <- as.data.frame(y0)
      y0 <- data.frame(cluster = rep(names(cluster_gene)[i], nrow(y0)), y0, stringsAsFactors = F)
      # filtering
      y <- subset(y0, Count >= 3 & p.adjust <= 0.05)
      return(y)
    })
    enrich_res <- do.call("rbind", enriched_LS)
    write.table(x = enrich_res, file = paste0(OUT, "/enrich_KEGG_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    res_in[[cell_type]][["enrich_KEGG"]] <- enrich_res
    # add stat
    res_in[[cell_type]]$cl_table_ftd$enrich_KEGG <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res$cluster
    res_in[[cell_type]]$cl_list_ftd$enrich_KEGG <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res$cluster
  }
  stopCluster(cl); rm(cl)
  return(res_in)
}

do_prepMSigDB <- function() {
  for(id in c("H", paste0("C", 1:7))) {
    cat(">", id, "\n")
    gs2gene <- msigdbr::msigdbr(species = "Homo sapiens", category = id)[, c("gs_name", "gene_symbol")]
    # build anno
    anno_build <- clusterProfiler:::build_Anno(path2gene = gs2gene, path2name = NA)
    saveRDS(object = anno_build, file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/MSigDB/MSigDB_", id, ".rds"))
  }
}

do_enrichMSigDB <- function(res_in, ctype = NULL, mdid = NULL, category = NULL, ncpu = 1) {
  # read gene sets
  if(is.null(category)) {
    category <- c("H", paste0("C", 1:7))[c(1:4,8)]
  }
  gs <- list()
  for(id in category) {
    gs[[id]] <- readRDS(file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/MSigDB/MSigDB_", id, ".rds"))
  }
  
  if(! is.null(ctype)) {
    cell_type_list <- ctype
  } else {
    cell_type_list <- names(res_in)
  }
  
  cl <- makeCluster(ncpu, type = "FORK")
  for(cell_type in cell_type_list) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    cluster_table_ftd <- res_in[[cell_type]]$cl_table_ftd
    cluster_table_ftd$cluster <- factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))
    cluster_gene <- split(cluster_table_ftd$gene, cluster_table_ftd$cluster)
    if(! is.null(mdid)) {
      cluster_gene <- cluster_gene[mdid]
    }
    enriched_LS <- parLapply(cl, seq_along(cluster_gene), function(i) {
      cat(">>", i, "\n")
      geneInput <- cluster_gene[[i]]
      enricher_LS <- lapply(gs, function(x_build) {
        y <- enricher_new(gene = geneInput, TERM2GENE = NA, TERM2NAME = NA, user_data = x_build, 
                          pAdjustMethod = "BH", minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 0.05)
        y <- as.data.frame(y)
        return(y)
      })
      allRes <- reshape2::melt(enricher_LS, measure.vars = NULL)
      allRes <- allRes[, c(ncol(allRes), 1:(ncol(allRes) - 1))]
      colnames(allRes)[1] <- "category"
      allRes <- data.frame(cluster = rep(names(cluster_gene)[i], nrow(allRes)), allRes, stringsAsFactors = F)
      # filtering
      allRes <- subset(allRes, Count >= 3 & p.adjust <= 0.05)
      return(allRes)
    })
    enrich_res <- do.call("rbind", enriched_LS)
    write.table(x = enrich_res, file = paste0(OUT, "/enrich_MSigDB_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    res_in[[cell_type]][["enrich_MSigDB"]] <- enrich_res
    # add stat
    res_in[[cell_type]]$cl_table_ftd$enrich_MSigDB <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res$cluster
    res_in[[cell_type]]$cl_list_ftd$enrich_MSigDB <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res$cluster
    # add stat (sub)
    enrich_res_TF <- subset(enrich_res, category == "C3" & (! grepl("_MIR", ID)))
    res_in[[cell_type]]$cl_table_ftd$enrich_TF <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res_TF$cluster
    res_in[[cell_type]]$cl_list_ftd$enrich_TF <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res_TF$cluster
    enrich_res_miRNA <- subset(enrich_res, category == "C3" & grepl("_MIR", ID))
    res_in[[cell_type]]$cl_table_ftd$enrich_miRNA <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res_miRNA$cluster
    res_in[[cell_type]]$cl_list_ftd$enrich_miRNA <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res_miRNA$cluster
  }
  stopCluster(cl); rm(cl)
  return(res_in)
}

do_prepTF <- function() {
  gs2gene <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/TF/TF_binding_rmdup.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1:2]
  dim(gs2gene)
  # gs_stat <- data.frame(num = table(gs2gene[, 1]), stringsAsFactors = F)
  # ggplot(gs_stat, aes(x = num.Freq)) + geom_histogram(fill = "lightblue", color = "blue", bins = 50) + xlab("Target gene number per TF") + ylab("TF number")
  # gs_stat <- data.frame(num = table(gs2gene[, 2]), stringsAsFactors = F)
  # ggplot(gs_stat, aes(x = num.Freq)) + geom_histogram(fill = "lightblue", color = "blue", bins = 50) + xlab("TF number per target gene") + ylab("Gene number")
  # build anno
  anno_build <- clusterProfiler:::build_Anno(path2gene = gs2gene, path2name = NA)
  saveRDS(object = anno_build, file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/TF/TF_", "binding", ".rds"))
}

do_prepMiRNA <- function() {
  gs2gene <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/miRNA/miRNA_binding_rmdup.txt", header = F, sep = "\t", stringsAsFactors = F)
  dim(gs2gene)
  # gs_stat <- data.frame(num = table(gs2gene[, 1]), stringsAsFactors = F)
  # ggplot(gs_stat, aes(x = num.Freq)) + geom_histogram(fill = "lightgreen", color = "darkgreen", bins = 50) + xlab("Target gene number per miRNA") + ylab("miRNA number")
  # gs_stat <- data.frame(num = table(gs2gene[, 2]), stringsAsFactors = F)
  # ggplot(gs_stat, aes(x = num.Freq)) + geom_histogram(fill = "lightgreen", color = "darkgreen", bins = 50) + xlab("miRNA number per target gene") + ylab("Gene number")
  # build anno
  anno_build <- clusterProfiler:::build_Anno(path2gene = gs2gene, path2name = NA)
  saveRDS(object = anno_build, file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/miRNA/miRNA_", "binding", ".rds"))
}

do_enrich <- function(res_in, ctype = NULL, mdid = NULL, db, category = NULL, expr_cutoff = NULL, res_in_sg = NULL, ncpu = 1) {
  # read gene sets
  if(is.null(category)) {
    category <- gsub("\\.rds", "", sapply(strsplit(list.files(path = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/", db), pattern = ".rds"), split = "_"), "[", 2))
  }
  gs <- list()
  for(id in category) {
    gs[[id]] <- readRDS(file = paste0("/rd/user/tianf/06-Human_cell_atlas/Data/", db, "/", db, "_", id, ".rds"))
  }
  
  if(! is.null(ctype)) {
    cell_type_list <- ctype
  } else {
    cell_type_list <- names(res_in)
  }
  
  cl <- makeCluster(ncpu, type = "FORK")
  for(cell_type in cell_type_list) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    cluster_table_ftd <- res_in[[cell_type]]$cl_table_ftd
    cluster_table_ftd$cluster <- factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))
    cluster_gene <- split(cluster_table_ftd$gene, cluster_table_ftd$cluster)
    if(! is.null(mdid)) {
      cluster_gene <- cluster_gene[mdid]
    }
    enriched_LS <- parLapply(cl, seq_along(cluster_gene), function(i) {
      cat(">>", i, "\n")
      geneInput <- cluster_gene[[i]]
      enricher_LS <- lapply(gs, function(x_build) {
        y <- enricher_new(gene = geneInput, TERM2GENE = NA, TERM2NAME = NA, user_data = x_build, 
                          pAdjustMethod = "BH", minGSSize = 1, maxGSSize = Inf, pvalueCutoff = 0.05)
        if(is.null(y)) {
          y <- as.data.frame(matrix(NA, nrow = 0, ncol = 9))
          colnames(y) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
        } else {
          y <- as.data.frame(y)
        }
        return(y)
      })
      allRes <- reshape2::melt(enricher_LS, measure.vars = NULL)
      allRes <- allRes[, c(ncol(allRes), 1:(ncol(allRes) - 1))]
      colnames(allRes)[1] <- "category"
      allRes <- data.frame(cluster = rep(names(cluster_gene)[i], nrow(allRes)), db = rep(db, nrow(allRes)), allRes, stringsAsFactors = F)
      # filtering
      allRes <- subset(allRes, Count >= 3 & p.adjust <= 0.05)
      if( (nrow(allRes) > 0) & (! is.null(expr_cutoff)) ) {
        if(cell_type != "merged") {
          allRes <- subset(allRes, ID %in% rownames(res_in[[cell_type]]$data))
          allRes$avgExpr <- rowMeans(res_in[[cell_type]]$data[allRes$ID, , drop = F])
          allRes <- subset(allRes, avgExpr >= expr_cutoff)
        } else {
          ct_high_avgCor <- colnames(res_in[[cell_type]]$avgCor)[res_in[[cell_type]]$avgCor[names(cluster_gene)[i], ] >= 0.088]
          res_in_sg <- res_in_sg[ct_high_avgCor]
          g_cand <- intersect(allRes$ID, rownames(res_in_sg[[1]]$data))
          if(length(g_cand) == 0) {
            allRes <- subset(allRes, F)
          } else if(length(g_cand) == 1) {
            expr_mean <- sapply(res_in_sg, function(x_md) { y <- rowMeans(x_md$data[g_cand, , drop = F]) })
            md_expressed <- any(expr_mean >= expr_cutoff)
            names(md_expressed) <- g_cand
            allRes <- subset(allRes, ID %in% names(md_expressed)[md_expressed])
          } else {
            expr_mean <- sapply(res_in_sg, function(x_md) { y <- rowMeans(x_md$data[g_cand, , drop = F]) })
            md_expressed <- apply(expr_mean, 1, function(x_expr) { any(x_expr >= expr_cutoff) })
            allRes <- subset(allRes, ID %in% names(md_expressed)[md_expressed])
          }
        }
      }
      return(allRes)
    })
    enrich_res <- do.call("rbind", enriched_LS)
    write.table(x = enrich_res, file = paste0(OUT, "/enrich_", db, "_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    res_in[[cell_type]][[paste0("enrich_", db)]] <- enrich_res
    # add stat
    res_in[[cell_type]]$cl_table_ftd[[paste0("enrich_", db)]] <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res$cluster
    res_in[[cell_type]]$cl_list_ftd[[paste0("enrich_", db)]] <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res$cluster
  }
  stopCluster(cl); rm(cl)
  return(res_in)
}

do_enrichr <- function(res_in, databases = c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")) {
  res_new <- res_in
  for(cell_type in names(res_in)) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    cluster_table_ftd <- res_in[[cell_type]]$cl_table_ftd
    cluster_gene <- split(cluster_table_ftd$gene, cluster_table_ftd$cluster)
    # no parallel
    enrich_res_LS <- lapply(cluster_gene, function(x) {
      y0 <- enrichr(genes = x, databases = databases)
      y <- reshape2:::melt(y0, measure.vars = NULL)
      colnames(y)[10] <- "database"
      return(y)
    })
    enrich_res <- reshape2:::melt(enrich_res_LS, measure.vars= NULL)
    colnames(enrich_res)[11] <- "cluster"
    enrich_res$Annotated <- as.numeric(gsub(".*/", "", enrich_res$Overlap))
    enrich_res$Overlap <- as.numeric(gsub("/.*", "", enrich_res$Overlap))
    enrich_res <- merge(enrich_res, unique(cluster_table_ftd[, c(1,3)]), by.x = "cluster", by.y = "cluster", sort = F)
    enrich_res <- merge(enrich_res, unique(cluster_table_ftd[, c(1,4)]), by.x = "cluster", by.y = "cluster", sort = F)
    # filtering
    enrich_res_ftd <- subset(enrich_res, Overlap>=3 & Adjusted.P.value <= 0.05)
    enrich_res_ftd <- enrich_res_ftd[, c(1,14,13,11,2,3,12,4:5,10)]
    write.table(x = enrich_res_ftd, file = paste0("03-expression/merged/geneModule/", samplingPos, "/geneEnrichment_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    res_new[[cell_type]][["enrichr_ftd"]] <- enrich_res_ftd
  }
  return(res_new)
}

do_plotEnrich <- function(res_in, ctype, mdid, db = "GO", category = NULL, shortName = T, charMaxLen = NULL, collapse = T, showNum = 10, upper.case = T, 
                          main = "auto", show_size = T, show.legend = F, bar.col = "#4DBBD5B2", asp.ratio = NULL, do.print = T, do.return = F) {
  enrich_res <- lapply(db, function(x_db) {
    y <- res_in[[ctype]][[paste0("enrich_", x_db)]]
    y <- y[, ! colnames(y) %in% c("db", "category")]
    y$db <- x_db
    return(y)
  })
  enrich_res <- do.call("rbind", enrich_res)
  if(is.null(enrich_res)) {
    return(NULL)
  }
  if(! is.null(category)) {
    env <- environment()
    enrich_res <- subset(enrich_res, category %in% get("category", env))
  }
  enrich_res <- subset(enrich_res, cluster == mdid, drop = F)
  enrich_res <- enrich_res[order(enrich_res$p.adjust), ]  # sort
  # reduce the char length
  if(shortName) {
    enrich_res$Description <- gsub(",.*", "", enrich_res$Description)
  }
  if(! is.null(charMaxLen)) {
    enrich_res$Description <- substr(enrich_res$Description, start = 1, stop = charMaxLen)
  }
  if(collapse) {
    enrich_res <- enrich_res[! duplicated(enrich_res$Description), ]
  }
  if(length(db) == 1) {
    enrich_res <- head(enrich_res, showNum)
  } else {
    enrich_res <- do.call("rbind", lapply(split(enrich_res, enrich_res$db)[db], function(x) { head(x, showNum) }))
  }
  if(upper.case) {
    enrich_res$Description <- Hmisc::capitalize(enrich_res$Description)
  }
  enrich_res$Description <- factor(enrich_res$Description, levels = rev(enrich_res$Description))
  
  if(! is.null(main)) {
    if(main == "auto") {
      main <- paste(ctype, mdid)
    }
    if(show_size) {
      md_size <- subset(res_in[[ctype]]$cl_list_ftd, cluster == mdid, "size", drop = T)
      main <- paste0(main, " (", md_size, "g)")
    }
  }
  
  if(nrow(enrich_res) == 0) {
    cat("No term will be plotted.\n")
  } else {
    gp <- ggplot(enrich_res, aes(x = Description, y = -log10(p.adjust), fill = db)) + geom_bar(stat = "identity", width = 0.9, show.legend = show.legend) + 
      coord_flip() + scale_fill_manual(values = bar.col) + 
      xlab(paste(paste(db, collapse = " / "), "terms")) + ylab(expression(paste(-log[10], " (FDR)"))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
      #scale_x_discrete(labels = rev(enrich_res$Description)) + 
      scale_y_continuous(expand = c(0, 0)) + 
      theme(axis.ticks.y = element_blank(), plot.margin = margin(7, 9, 7, 7, "pt")) + ggtitle(main)
    if(length(db) > 1) {
      gp <- gp + facet_grid(db ~ ., scales = "free_y", space = "free_y") + theme(strip.background = element_blank())
    }
    if(nrow(enrich_res) < 5) {
      gp <- gp + theme(aspect.ratio = 0.2)
    }
    if(! is.null(asp.ratio)) {
      gp <- gp + theme(aspect.ratio = asp.ratio)
    }
    if(do.print) {
      print(gp)
    } else if(do.return) {
      return(gp)
    }
  }
}

do_plotModuleMap <- function(res_in, ctype, min.avgExpr = -Inf, min.avgCor = -Inf, sortBy = "avgCor", mdNum = NULL, 
                             showModuleID = T, widths = c(0.9, 0.65, 0.5, 0.6, 1)) {
  cl_list_ftd <- res_in[[ctype]]$cl_list_ftd
  # filtering
  cl_list_ftd <- subset(cl_list_ftd, avgExpr >= min.avgExpr & avgCor >= min.avgCor)
  # sort
  cl_list_ftd <- cl_list_ftd[order(cl_list_ftd[[sortBy]], decreasing = T), ]
  if(! is.null(mdNum)) {
    cl_list_ftd <- head(cl_list_ftd, mdNum)
  }
  cl_list_ftd$cluster <- factor(cl_list_ftd$cluster, levels = rev(cl_list_ftd$cluster))
  # gene stat
  cl_list_ftd_gstat_melted <- reshape2::melt(cl_list_ftd[, c("cluster", "protein_coding", "lncRNA", "sncRNA", "pseudogene", "Ig/TcR")], id.vars = "cluster")
  levels(cl_list_ftd_gstat_melted$variable) <- c("Pc", "Lnc", "Snc", "Pse", "Ig/TR")
  # anno stat
  cl_list_ftd_melted <- reshape2::melt(cl_list_ftd[, c("cluster", "enrich_TF", "enrich_miRNA", "enrich_GO", "enrich_KEGG", "enrich_PPI")], id.vars = "cluster")
  cl_list_ftd_melted$variable <- gsub("^enrich_", "", cl_list_ftd_melted$variable)
  cl_list_ftd_melted$variable <- factor(cl_list_ftd_melted$variable, levels = unique(cl_list_ftd_melted$variable))

  p1 <- ggplot(cl_list_ftd, aes(x = 1, y = cluster)) + geom_tile(aes(fill = avgCor)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
          legend.direction = "horizontal", legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = -10), 
          aspect.ratio = 5) + 
    labs(fill = "Correlation") + ylab(paste0("Module (N=", nrow(cl_list_ftd), ")")) + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
    scale_x_continuous(expand = c(0, 0)) + scale_fill_gradientn(colours = colorRampPalette(c("white", "red"))(100), limits = c(0, 0.3), breaks = c(0, 0.15, 0.3))
  if(! showModuleID) {
    p1 <- p1 + theme(axis.text.y = element_blank())
  }
  
  p2 <- ggplot(cl_list_ftd, aes(x = 1, y = cluster)) + geom_tile(aes(fill = log2(avgExpr + 1))) + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
          legend.direction = "horizontal", legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = -10), legend.title = element_text(margin = margin(b = -2)), 
          aspect.ratio = 5) + 
    labs(fill = expression(log[2] ~ "(" * Expr * "+1)")) + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
    scale_x_continuous(expand = c(0, 0)) + scale_fill_gradientn(colours = colorRampPalette(c("white","palegreen3"))(100))
  
  p3 <- ggplot(cl_list_ftd, aes(x = size, y = cluster)) + geom_segment(aes(x = 0, xend = size, y = cluster, yend = cluster)) + geom_point(color = "grey60", size = 1) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          plot.margin = margin(t = 8, b = 28)) + 
    xlab("Module size") + coord_cartesian(clip = "off") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(cl_list_ftd$size) * 1.02))
  
  p4 <- ggplot(cl_list_ftd_gstat_melted, aes(x = cluster, y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = position_stack(reverse = T)) + geom_bar(data = cl_list_ftd, aes(y = tf), fill = "brown1", stat = "identity") + 
    coord_flip(clip = "off") + 
    theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          legend.direction = "horizontal", legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = -12, b = -2)) + 
    guides(fill = guide_legend(ncol = 2)) + 
    labs(fill = NULL) + 
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) + scale_fill_manual(values = c("hotpink", "dodgerblue", "lightskyblue", "grey60", "grey90"))
    
  p5 <- ggplot(cl_list_ftd_melted, aes(x = 1, y = cluster)) + geom_tile(aes(fill = value)) + facet_grid(. ~ variable, switch = "x") + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
          legend.direction = "horizontal", legend.position = "bottom", legend.justification = "center", legend.margin = margin(t = -10, b = 18), 
          strip.background = element_blank(), strip.text.x = element_text(margin = margin(t = 10))) + 
    labs(fill = NULL) + 
    scale_x_continuous(expand = c(0, 0)) + scale_fill_manual(values = c("grey80", "skyblue"))
  p5 <- ggplotGrob(p5)
  for(i in which(grepl("strip-b", p5$layout$name))){ p5$grobs[[i]]$layout$clip <- "off" }
  
  grid.arrange(p1, p2, p3, p4, p5, nrow = 1, widths = widths, top = ctype)
}

do_prepPPI <- function() {
  # create string_db
  library("STRINGdb")
  dbDir <- "/rd/user/tianf/06-Human_cell_atlas/Data/STRINGdb"
  string_db <- STRINGdb$new(version = "10", species = 9606, score_threshold = 400, input_directory = dbDir)
  string_db$load()
  saveRDS(string_db, file = paste0(dbDir, "/string_db.rds"))
  
  # create ID mapping file for all genes
  id_DF <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/gene_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)
  colnames(id_DF) <- c("ensembl_id", "gene")
  # method1
  id_mapped1 <- string_db$map(id_DF, "gene", removeUnmappedRows = T)
  write.table(x = id_mapped1, file = paste0(dbDir, "/id_mapping_pkg.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  # method2
  # protein_aliases <- read.table(paste0(dbDir, "/id_mapping.txt"), header = F, sep = "\t", stringsAsFactors = F)
  # dim(protein_aliases)
  # colnames(protein_aliases) <- c("STRING_id", "alias", "source")
  # id_mapped2 <- merge(id_DF, protein_aliases, by.x = "gene", by.y = "alias", sort = F, all.x = F)[, -6]
  # # cmp
  # length(setdiff(id_mapped1$gene, id_mapped2$gene))
  # id_mappeds <- merge(id_mapped1, id_mapped2, by = "gene", sort = F)
  # all(id_mappeds$STRING_id.x==id_mappeds$STRING_id.y)
}

do_enrichPPI <- function(res_in, ctype = NULL, mdid = NULL, ncpu = 1) {
  # get string_db
  string_db <- readRDS(file = "/rd/user/tianf/06-Human_cell_atlas/Data/STRINGdb/string_db.rds")
  # get ID mapping info
  gene2id <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/STRINGdb/id_mapping_pkg.txt", header = T, sep = "\t", stringsAsFactors = F)
  
  if(! is.null(ctype)) {
    cell_type_list <- ctype
  } else {
    cell_type_list <- names(res_in)
  }
  
  cl <- makeCluster(ncpu, type = "FORK")
  for(cell_type in cell_type_list) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    cluster_table_ftd <- res_in[[cell_type]]$cl_table_ftd
    cluster_table_ftd$cluster <- factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))
    cluster_gene <- split(cluster_table_ftd$gene, cluster_table_ftd$cluster)
    if(! is.null(mdid)) {
      cluster_gene <- cluster_gene[mdid]
    }
    enriched_LS <- parLapply(cl, seq_along(cluster_gene), function(i) {
      cat(">>", i, "\n")
      geneInput <- cluster_gene[[i]]
      # id mapping
      id_DF <- data.frame(gene = geneInput, stringsAsFactors = F)
      id_mapped <- merge(id_DF, gene2id, by = "gene", sort = F)
      id_mapped_mulMp <- id_mapped$gene[duplicated(id_mapped$gene)]
      num_mt <- nrow(id_DF)
      num_m0 <- nrow(id_DF) - length(unique(id_mapped$gene))
      num_m1 <- length(unique(id_mapped$gene)) - length(unique(id_mapped_mulMp))
      num_mm <- length(id_mapped_mulMp)
      rat_m1 <- num_m1 / num_mt
      # rmdup
      id_mapped <- id_mapped[! duplicated(id_mapped$STRING_id), ]
      num_pro <- length(id_mapped$STRING_id)
      # calc enrichment
      hitsWithEdges <- id_mapped$STRING_id[id_mapped$STRING_id %in% igraph::V(string_db$graph)$name]
      enrichment <- string_db$get_ppi_enrichment(id_mapped$STRING_id)
      num_exp <- enrichment$lambda
      num_obs <- length(igraph::E(igraph::induced.subgraph(string_db$graph, hitsWithEdges)))
      pvalue <- enrichment$enrichment
      
      # if(do_plot) {
      #   #string_db$plot_network(id_mapped$STRING_id, add_summary = F, add_link = F)
      #   img <- string_db$get_png(id_mapped$STRING_id)
      #   par(mar = c(0, 0, 0, 0))
      #   plot(1:(dim(img)[2]), type = "n", xaxt = "n", yaxt = "n", 
      #        xlab = "", ylab = "", ylim = c(1, dim(img)[1]), xlim = c(1, (dim(img)[2])), asp = 1, bty='n')
      #   rasterImage(img, 1, 1, dim(img)[2], dim(img)[1])
      # }
      # if(get_link) {
      #   link_text <- string_db$get_link(id_mapped$STRING_id)
      #   print(link_text)
      # }
      
      y0 <- data.frame(total = num_mt, no_hit = num_m0, single_hit = num_m1, 
                       num_pro = num_pro, num_exp = num_exp, num_obs = num_obs, pvalue = pvalue, 
                       geneID = paste(geneInput, collapse = "/"), STRING_id = paste(id_mapped$STRING_id, collapse = "/"), stringsAsFactors = F)
      y0 <- data.frame(cluster = rep(names(cluster_gene)[i], nrow(y0)), y0, stringsAsFactors = F)
      # filering
      y <- subset(y0, pvalue <= 0.05)
      return(y)
    })
    enrich_res <- do.call("rbind", enriched_LS)
    write.table(x = enrich_res, file = paste0(OUT, "/enrich_PPI_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    res_in[[cell_type]][["enrich_PPI"]] <- enrich_res
    # add stat
    res_in[[cell_type]]$cl_table_ftd$enrich_PPI <- res_in[[cell_type]]$cl_table_ftd$cluster %in% enrich_res$cluster
    res_in[[cell_type]]$cl_list_ftd$enrich_PPI <- res_in[[cell_type]]$cl_list_ftd$cluster %in% enrich_res$cluster
  }
  stopCluster(cl); rm(cl)
  return(res_in)
}

do_mergeModule <- function(res_in, ov_cutoff = 0.9, rename = T, verbose = F) {
  # read gene type info
  gene_type <- read.table("/rd/user/tianf/06-Human_cell_atlas/Genomes/human/gene_type_class.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 2)
  colnames(gene_type) <- c("ensembl_id", "type", "class")
  gene_type$class <- factor(gene_type$class, levels = unique(gene_type$class)[c(4,3,2,1,5)])
  # read tf list
  tf_list <- read.table("/rd/user/tianf/06-Human_cell_atlas/Data/AnimalTFDB/human_TF_list.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]
  
  # joint module analysis
  cl_LS <- lapply(res_in, function(x) {
    cl_tb <- x$cl_table_ftd
    cl_tb$cluster <- factor(cl_tb$cluster, levels = unique(cl_tb$cluster))
    cl_sl <- split(x = cl_tb$gene, f = cl_tb$cluster)
    return(cl_sl)
  })
  
  cl_rmdup <- list()
  for(i in names(cl_LS)) {  # for each cell type
    cat("> ", i, "\n")
    if(length(cl_rmdup) == 0) {
      cl_rmdup <- cl_LS[[i]]
      names(cl_rmdup) <- paste0(i, "_", names(cl_rmdup))
      next
    }
    
    #cl_melted <- melt(cl_LS[[i]])
    #colnames(cl_melted) <- c("gene", "module")
    
    for(j in names(cl_LS[[i]])) { # for each module in this cell type (from 2nd)
      kc <- 0
      g_query <- cl_LS[[i]][[j]]
      for(k in names(cl_rmdup)) { # for each md in rmdup
        g_ref <- cl_rmdup[[k]]
        gene_it <- intersect(g_ref, g_query)
        gene_un <- union(g_ref, g_query)
        num_ia <- length(g_query)
        num_ib <- length(g_ref)
        num_it <- length(gene_it)
        num_un <- length(gene_un)
        rat_ia <- num_it / num_ia
        rat_ib <- num_it / num_ib
        rat_it <- num_it / num_un
        if( rat_it >= ov_cutoff ) {
          if(verbose) {
            cat(i, j, k, "|", num_ia, num_ib, num_it, num_un, "|", rat_ia, rat_ib, rat_it, "\n")
          }
          cl_rmdup[[k]] <- gene_it
          kc <- 1
        }
      }
      if(kc == 0) {
        newid <- paste0(i, "_", j)
        cl_rmdup[[newid]] <- g_query
      }
    }
  }
  if(rename) {
    names(cl_rmdup) <- paste0("MD", seq_along(cl_rmdup))
  }
  
  cluster_table <- melt(cl_rmdup)
  colnames(cluster_table) <- c("gene", "cluster")
  cluster_table$gene <- as.character(cluster_table$gene)
  cluster_gene <- split(cluster_table$gene, cluster_table$cluster)[unique(cluster_table$cluster)]
  # add size
  cluster_size <- data.frame(size = table(cluster_table$cluster), stringsAsFactors = F)
  colnames(cluster_size) <- c("cluster", "size")
  cluster_table <- merge(cluster_table, cluster_size, by.x = 2, by.y = 1, sort = F)
  # add boundary
  cluster_bdl <- c(1, cumsum(cluster_table[! duplicated(cluster_table$cluster), "size"])[-length(unique(cluster_table$cluster))] + 1)
  cluster_bdr <- cumsum(cluster_table[! duplicated(cluster_table$cluster), "size"])
  cluster_bd <- data.frame(cluster=unique(cluster_table$cluster), boundary=paste(cluster_bdl, cluster_bdr, sep = ":"), stringsAsFactors = F)
  cluster_table <- merge(cluster_table, cluster_bd, by = "cluster", sort = F)
  # add expr stat
  cluster_exprStat <- t(sapply(cluster_gene, function(x) { ot <- sapply(res_in, function(k) { y0 <- k$data[x, ]; y <- mean(unlist(y0), na.rm = T); return(y) }); return(ot) }))
  #cluster_table <- merge(cluster_table, cluster_exprStat, by.x = "cluster", by.y = 0, sort = F)
  # add cor stat
  cluster_corStat <- t(sapply(cluster_gene, function(x) { ot <- sapply(res_in, function(k) { xs <- intersect(x, rownames(k$cor_cld)); y0 <- k$cor_cld[xs, xs]; y <- mean(y0[upper.tri(y0)]); return(y) }); return(ot) }))
  #cluster_table <- merge(cluster_table, cluster_corStat, by.x = "cluster", by.y = 0, sort = F)
  # add gene type stat
  cluster_gtStat <- t(sapply(cluster_gene, function(x) { y <- table(gene_type[x, "class"]) / length(x); return(y) }))
  cluster_table <- merge(cluster_table, cluster_gtStat, by.x = "cluster", by.y = 0, sort = F)
  # add tf stat
  cluster_tfStat <- data.frame(tf = sapply(cluster_gene, function(x) { y <- sum(x %in% tf_list) / length(x); return(y) }))
  cluster_table <- merge(cluster_table, cluster_tfStat, by.x = "cluster", by.y = 0, sort = F)
  # filtering
  cluster_table_ftd <- cluster_table
  
  mes_in <- list()
  mes_in[["merged"]][["cl_table_ftd"]] <- cluster_table_ftd
  cluster_list_ftd <- unique(cluster_table_ftd[, -2])
  rownames(cluster_list_ftd) <- NULL
  mes_in[["merged"]][["cl_list_ftd"]] <- cluster_list_ftd
  mes_in[["merged"]][["avgExpr"]] <- cluster_exprStat
  mes_in[["merged"]][["avgCor"]] <- cluster_corStat
  return(mes_in)
}
