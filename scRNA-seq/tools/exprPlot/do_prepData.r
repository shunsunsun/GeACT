
setwd("~/lustre/06-Human_cell_atlas/tools/exprPlot/")

suppressMessages(library("arrow"))
library("parallel")

OUT <- "data"
dir.create(OUT, showWarnings = F, recursive = T)

# load gene expression matrix ( log2(CPM+1) )
exprMatr <- read_feather(file = "~/lustre/06-Human_cell_atlas/pooled_data_all/All/03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.feather")
exprMatr <- as.data.frame(exprMatr)
rownames(exprMatr) <- read.table(file = "~/lustre/06-Human_cell_atlas/pooled_data_all/All/03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.gene", header = F, sep = "\t", stringsAsFactors = F)$V1

# load cell meta
cellMeta <- read.table("~/lustre/06-Human_cell_atlas/pooled_data_all/All/cell_metatable_filtered_aligned.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
all(colnames(exprMatr) == rownames(cellMeta))

ts_all <- unique(cellMeta$tissue)

cl <- makeCluster(10, type = "FORK")
parLapply(cl, ts_all, function(x_ts) {
  cat(">", x_ts, "\n")
  cellMeta_sub <- subset(cellMeta, tissue == x_ts, c("ident", "group", "stage", "tSNE_1_ali", "tSNE_2_ali", "UMAP_1_ali", "UMAP_2_ali"))
  exprMatr_sub <- exprMatr[, rownames(cellMeta_sub)]
  cellMatr <- cbind(cellMeta_sub, t(exprMatr_sub))
  write_feather(x = cellMatr, sink = paste0(OUT, "/", gsub(" ", "_", x_ts), "_cellMatr.feather"))
  write.table(x = rownames(cellMatr), file = paste0(OUT, "/", gsub(" ", "_", x_ts), "_cellMatr.cell"), row.names = F, col.names = F, quote = F, sep = "\t")
  #write_feather(x = exprMatr_sub, sink = paste0(OUT, "/", gsub(" ", "_", x_ts), "_exprMatr.feather"))
  #write.table(x = rownames(exprMatr_sub), file = paste0(OUT, "/", gsub(" ", "_", x_ts), "_exprMatr.gene"), row.names = F, col.names = F, quote = F, sep = "\t")
  #write.table(x = cellMeta_sub, file = paste0(OUT, "/", gsub(" ", "_", x_ts), "_cellMeta.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
})
stopCluster(cl); rm(cl)
