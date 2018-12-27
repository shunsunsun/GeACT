md.preprocess <- function(expr_data, cellMetaData, outdir = "03-expression/merged/modules") {
  # cell filtering and normalization
  dir.create(path = outdir, showWarnings = F, recursive = T)
  # CPM
  expr_data <- sweep(x = expr_data, MARGIN = 2, STATS = colSums(expr_data), FUN = "/") * 0.1 * 1e6 # scale by 0.1M
  # filtering genes
  cat("Genes before filtering:", nrow(expr_data), "\n")
  expr_data <- expr_data[rowSums(expr_data>0)/ncol(expr_data)>0.05, ]
  cat("Genes after filtering:", nrow(expr_data), "\n")
  # split the gene expression table by cell type
  cellCluster_LS <- split(rownames(cellMetaData), cellMetaData$expr.ident)
  lengths(cellCluster_LS)

  expr_sub_LS <- lapply(cellCluster_LS, function(x) {
    y <- expr_data[, colnames(expr_data)%in%x]
    return(y)
  })
  return(expr_sub_LS)
}

md.clustering <- function(expr_sub_LS) {
  # HC clustering
  hc_LS <- list()
  cor_LS <- list()
  for(cell_type in names(expr_sub_LS)) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    expr_input <- expr_sub_LS[[cell_type]]
    # HC clustering
    expr_cor <- cor(t(expr_input))
    ds <- as.dist(1-expr_cor)
    hc <- hclust(ds, method = "average")
    expr_cor_clustered <- expr_cor[hc$order, hc$order]
    hc_LS[[cell_type]] <- hc
    cor_LS[[cell_type]] <- expr_cor_clustered
    rm(cell_type, cell_type_label, expr_input, expr_cor, ds, hc, expr_cor_clustered)
  }
  return(list(hc=hc_LS, cor=cor_LS))
}

md.module <- function(hc_cor_LS, cor_cutoff = 0.225, outdir = "03-expression/merged/modules") {
  # Gene module detection
  # cat(">> Gene module detection\n")
  cluster_table_LS <- list()
  hc_LS <- hc_cor_LS[["hc"]]
  cor_LS <- hc_cor_LS[["cor"]]
  for(cell_type in names(hc_LS)) {
    cat(">", cell_type, "\n")
    cell_type_label <- gsub(" ", "_", cell_type)
    hc <- hc_LS[[cell_type]]
    expr_cor_clustered <- cor_LS[[cell_type]]
    # cut tree
    cluster_table <- data.frame(gene = names(cutree(tree = hc, h = 1-cor_cutoff)), cluster = cutree(tree = hc, h = 1-cor_cutoff), stringsAsFactors = F)
    cluster_table$cluster <- paste0("M", cluster_table$cluster)
    # add size
    cluster_size <- data.frame(size = table(cluster_table$cluster))
    colnames(cluster_size) <- c("cluster", "size")
    cluster_table <- merge(cluster_table, cluster_size, by.x = 2, by.y = 1, sort = F)
    # add boundary
    cluster_table <- merge(data.frame(gene=rownames(expr_cor_clustered), stringsAsFactors = F), cluster_table, by.x = 1, by.y = 2, sort = F) # sort by gene ID in clustered
    cluster_bdl <- c(1, cumsum(cluster_table[! duplicated(cluster_table$cluster), "size"])[-length(unique(cluster_table$cluster))] + 1)
    cluster_bdr <- cumsum(cluster_table[! duplicated(cluster_table$cluster), "size"])
    cluster_bd <- data.frame(cluster=unique(cluster_table$cluster), boundary=paste(cluster_bdl, cluster_bdr, sep = "-"))
    cluster_table <- merge(cluster_table, cluster_bd, by.x = "cluster", by.y = "cluster", sort = F)
    # filtering by size
    cat("Raw:", length(unique(cluster_table$cluster)), "modules with", nrow(cluster_table), "genes\n")
    cluster_table_ftd <- subset(cluster_table, size>=10)
    cat("Clean:", length(unique(cluster_table_ftd$cluster)), "modules with", nrow(cluster_table_ftd), "genes\n")
    write.table(x = cluster_table_ftd, file = paste0(outdir, "/geneCluster_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    cluster_table_LS[[cell_type]] <- cluster_table_ftd
    # filtering correlation matrix
    expr_cor_ftd <- expr_cor_clustered[rownames(expr_cor_clustered)%in%cluster_table_ftd$gene, colnames(expr_cor_clustered)%in%cluster_table_ftd$gene]

    png(paste0(outdir, "/geneCorrelation_", cell_type_label, ".png"), height = 1600, width = 1600, res = 300)
    # all modules
    ph <- pheatmap(mat = expr_cor_ftd, cluster_rows = F, cluster_cols = F, breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1),
             color = c("blue",colorRampPalette(c("blue", "white", "red"))(100),"red"), show_rownames = F, show_colnames = F, legend = T,
             main = cell_type, silent = T)
    print(ph)
    dev.off()

    rm(cell_type, cell_type_label, hc, expr_cor_clustered,
       cor_cutoff, cluster_table, cluster_table_ftd, cluster_size, cluster_bdl, cluster_bdr, cluster_bd,
       expr_cor_ftd)
  }
  return(cluster_table_LS)
}
