# gene module detection

# demo START >>>>>>>>>>>>>>>>>>>>>>
library("pheatmap")

# Load the dataset
#detach("package:gi", unload=TRUE)
library("gi")

projectID <- "geact2"
sessionID <- "5c2325197e71c30a28f5359b"

expr_data <- t(getExprMatrix(projectID, sessionID))
dim(expr_data)
cellStat <- getMeta(projectID, sessionID)
dim(cellStat)

# CPM
expr_data <- sweep(x = expr_data, MARGIN = 2, STATS = colSums(expr_data), FUN = "/") * 0.1 * 1e6 # scale by 0.1M

# only show mono cells
expr_data <- expr_data[, grep("^Mono", colnames(expr_data), value = T)]
expr_data <- expr_data[rowSums(expr_data>0)/ncol(expr_data)>0.05, ]
dim(expr_data)
cellMetaData <- data.frame(cellStat[colnames(expr_data), ], expr.ident = "Mono", stringsAsFactors = F)
dim(cellMetaData)

# split the gene expression table by cell type
cellCluster_LS <- split(rownames(cellMetaData), cellMetaData$expr.ident)
lengths(cellCluster_LS)

expr_sub_LS <- lapply(cellCluster_LS, function(x) {
  y <- expr_data[, colnames(expr_data)%in%x]
  return(y)
})
names(expr_sub_LS)

dir.create(path = "03-expression/merged/modules", showWarnings = F, recursive = T)

# 1. HC clustering ---
cat(">> HC clustering\n")
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

# 2. Gene module detection ----
cat(">> Gene module detection\n")
cluster_table_LS <- list()
for(cell_type in names(expr_sub_LS)) {
  cat(">", cell_type, "\n")
  cell_type_label <- gsub(" ", "_", cell_type)
  hc <- hc_LS[[cell_type]]
  expr_cor_clustered <- cor_LS[[cell_type]]
  # cut tree
  cor_cutoff <- 0.225  # XXX
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
  write.table(x = cluster_table_ftd, file = paste0("03-expression/merged/modules/geneCluster_", cell_type_label, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  cluster_table_LS[[cell_type]] <- cluster_table_ftd
  # filtering correlation matrix
  expr_cor_ftd <- expr_cor_clustered[rownames(expr_cor_clustered)%in%cluster_table_ftd$gene, colnames(expr_cor_clustered)%in%cluster_table_ftd$gene]
  #write.table(x = expr_cor_ftd, file = paste0("03-expression/merged/modules/geneCorrelation_", cell_type_label, ".txt", row.names = T, col.names = NA, quote = F, sep = "\t")
  
  png(paste0("03-expression/merged/modules/geneCorrelation_", cell_type_label, ".png"), height = 1600, width = 1600, res = 300)
  # all modules
  pheatmap(mat = expr_cor_ftd, cluster_rows = F, cluster_cols = F, breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1),
           #gaps_row = cumsum(cluster_table_ftd[! duplicated(cluster_table_ftd$cluster), "size"])[1],
           #gaps_col = cumsum(cluster_table_ftd[! duplicated(cluster_table_ftd$cluster), "size"])[1],
           color = c("blue",colorRampPalette(c("blue", "white", "red"))(100),"red"), show_rownames = F, show_colnames = F, legend = T,
           main = cell_type)
  dev.off()
  
  # pdf(paste0("03-expression/merged/modules/geneCorrelation_", cell_type_label, ".pdf"), height = 8, width = 8)
  # # relationship between modules
  # expr_cor_rowMeans <- aggregate(x = expr_cor_ftd, by = list(factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))), mean)
  # rownames(expr_cor_rowMeans) <- expr_cor_rowMeans[, 1]
  # expr_cor_rowMeans <- expr_cor_rowMeans[, -1]
  # expr_cor_mean <- aggregate(x = t(expr_cor_rowMeans), by = list(factor(cluster_table_ftd$cluster, levels = unique(cluster_table_ftd$cluster))), mean)
  # rownames(expr_cor_mean) <- expr_cor_mean[, 1]
  # expr_cor_mean <- expr_cor_mean[, -1]
  # 
  # display_numbers_MT <- matrix(ifelse(abs(expr_cor_mean) > cor_cutoff, "*", ""), nrow(expr_cor_mean))
  # diag(display_numbers_MT) <- ""
  # pheatmap(mat = expr_cor_mean, cluster_rows = F, cluster_cols = F, breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1), 
  #          color = c("blue",colorRampPalette(c("blue", "white", "red"))(100),"red"), show_rownames = T, show_colnames = T, legend = F, 
  #          display_numbers = display_numbers_MT, fontsize = 5, number_color = "black", main = cell_type)
  # 
  # # the most negative values
  # expr_cor_mean_melted <- melt(as.matrix(expr_cor_mean))
  # expr_cor_mean_melted$id1 <- as.numeric(gsub("^M", "", expr_cor_mean_melted$Var1))
  # expr_cor_mean_melted$id2 <- as.numeric(gsub("^M", "", expr_cor_mean_melted$Var2))
  # subset(expr_cor_mean_melted, value<(-cor_cutoff) & id1<id2)
  # 
  # do_cluster2bd <- function(x) {
  #   y0 <- subset(cluster_bd, cluster==x, select = "boundary", drop = T)
  #   #print(y0)
  #   ya <- as.numeric(gsub("-.*", "", y0))
  #   yb <- as.numeric(gsub(".*-", "", y0))
  #   y <- ya:yb
  #   return(y)
  # }
  # 
  # query_cluster <- as.character(unique(unlist(subset(expr_cor_mean_melted, value<(-cor_cutoff) & id1<id2)[, 1:2])))
  # cat("Modules with negative cor:", query_cluster, "\n")
  # 
  # if(length(query_cluster)>0) {
  #   subid <- c()
  #   subsize <- c()
  #   for(i in query_cluster) {
  #     subid <- c(subid, do_cluster2bd(i))
  #     subsize <- c(subsize, subset(cluster_size, cluster==i, select = "size", drop = T))
  #   }
  #   annotation_row_DF <- data.frame(row.names = rownames(expr_cor_clustered)[subid], module = cluster_table$cluster[subid])
  #   annotation_row_DF$module <- factor(annotation_row_DF$module, levels = unique(annotation_row_DF$module))
  #   pheatmap(mat = expr_cor_clustered[subid, subid], cluster_rows = F, cluster_cols = F, breaks = c(-1, seq(from = -0.6, to = 0.6, length.out = 100), 1), 
  #            gaps_row = cumsum(subsize), 
  #            gaps_col = cumsum(subsize), 
  #            color = c("blue",colorRampPalette(c("blue", "white", "red"))(100),"red"), show_rownames = F, show_colnames = F, legend = F, 
  #            annotation_row = annotation_row_DF, annotation_legend = T, annotation_names_row = F, 
  #            main = cell_type)
  #   rm(subid, subsize, i, annotation_row_DF)
  # } else {
  #   #cat("No negative inter-module correlations\n")
  # }
  # dev.off()
  rm(cell_type, cell_type_label, hc, expr_cor_clustered, 
     cor_cutoff, cluster_table, cluster_table_ftd, cluster_size, cluster_bdl, cluster_bdr, cluster_bd, 
     expr_cor_ftd)
  # rm(expr_cor_rowMeans, expr_cor_mean, expr_cor_mean_melted, 
  #    do_cluster2bd, display_numbers_MT, query_cluster)
}

