# pool datasets from different seqID
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

library("parallel")
suppressMessages(library("arrow"))

# options
dts <- data.frame(tissue = list.files(path = "..", pattern = "^[0-9]"), stringsAsFactors = F)
#out <- "X"

# 1.1 UMI (before filtering)
cl <- makeCluster(min(nrow(dts), 10), type = "FORK")
dts_LS <- parLapply(cl, split(dts, 1:nrow(dts)), function(x) { 
  print(x)
  y <- read_feather(file.path("..", x[1], "03-expression/merged/filtering/UMIcount_unfiltered.feather"))
  y <- as.data.frame(y)
  y_gene <- read.table(file = file.path("..", x[1], "03-expression/merged/filtering/UMIcount_unfiltered.gene"), header = F, sep = "\t", stringsAsFactors = F)
  rownames(y) <- y_gene$V1
  return(y)
})
stopCluster(cl); rm(cl)
names(dts_LS) <- NULL
dts_DF <- do.call("cbind", dts_LS)
dir.create(path = paste0("03-expression/merged/filtering"), showWarnings = F, recursive = T)
# write.table vs fwrite (1000 rows * 31392 columns: 186.282s ~ 4.647s)
data.table::fwrite(x = dts_DF, file = paste0("03-expression/merged/filtering/UMIcount_unfiltered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_unfiltered.txt", " > ", "03-expression/merged/filtering/UMIcount_unfiltered.txt.gz"))
# for quick read
write_feather(x = dts_DF, sink = paste0("03-expression/merged/filtering/UMIcount_unfiltered.feather"))
write.table(x = rownames(dts_DF), file = paste0("03-expression/merged/filtering/UMIcount_unfiltered.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# 1.2 UMI (after filtering)
# filtering cells
ftc_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/03-expression/merged/filtering/filtering_cells.txt"), header = T, sep = "\t", stringsAsFactors = F)
  return(y)
})
names(ftc_LS) <- NULL
ftc_DF <- do.call("rbind", ftc_LS)
colnames(ftc_DF)[1] <- "cell"
all(colnames(dts_DF) == ftc_DF$cell)
dir.create(path = paste0("03-expression/merged/filtering"), showWarnings = F, recursive = T)
data.table::fwrite(x = ftc_DF, file = paste0("03-expression/merged/filtering/filtering_cells.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
dts_cellftd <- dts_DF[, subset(ftc_DF, filter, "cell", drop = T)]

data.table::fwrite(x = dts_cellftd, file = paste0("03-expression/merged/filtering/UMIcount_cellFiltered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_cellFiltered.txt", " > ", "03-expression/merged/filtering/UMIcount_cellFiltered.txt.gz"))
# for quick read
write_feather(x = dts_cellftd, sink = paste0("03-expression/merged/filtering/UMIcount_cellFiltered.feather"))
write.table(x = rownames(dts_cellftd), file = paste0("03-expression/merged/filtering/UMIcount_cellFiltered.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# norm (using the genes before filtering)
dts_cellftd_CPM <- as.data.frame(apply(dts_cellftd, 2, function(x) (x / sum(x))) * 1e6)
data.table::fwrite(x = dts_cellftd_CPM, file = paste0("03-expression/merged/filtering/UMIcount_cellFiltered_CPM.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_cellFiltered_CPM.txt", " > ", "03-expression/merged/filtering/UMIcount_cellFiltered_CPM.txt.gz"))
# for quick read
write_feather(x = dts_cellftd_CPM, sink = paste0("03-expression/merged/filtering/UMIcount_cellFiltered_CPM.feather"))
write.table(x = rownames(dts_cellftd_CPM), file = paste0("03-expression/merged/filtering/UMIcount_cellFiltered_CPM.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# log2 norm
dts_cellftd_log2CPM <- log2(dts_cellftd_CPM + 1)
data.table::fwrite(x = dts_cellftd_log2CPM, file = paste0("03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.txt", " > ", "03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.txt.gz"))
# for quick read
write_feather(x = dts_cellftd_log2CPM, sink = paste0("03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.feather"))
write.table(x = rownames(dts_cellftd_log2CPM), file = paste0("03-expression/merged/filtering/UMIcount_cellFiltered_log2CPM.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# filtering genes
nCell_expressed <- rowSums(dts_cellftd > 0)
dts_ftd <- dts_cellftd[nCell_expressed >= 10, ]
dim(dts_ftd)
data.table::fwrite(x = dts_ftd, file = paste0("03-expression/merged/filtering/UMIcount_filtered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_filtered.txt", " > ", "03-expression/merged/filtering/UMIcount_filtered.txt.gz"))
# for quick read
write_feather(x = dts_ftd, sink = paste0("03-expression/merged/filtering/UMIcount_filtered.feather"))
write.table(x = rownames(dts_ftd), file = paste0("03-expression/merged/filtering/UMIcount_filtered.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# norm (using the genes after filtering)
dts_ftd_CPM <- sweep(x = dts_ftd, MARGIN = 2, STATS = colSums(dts_ftd), FUN = "/") * 1e6
data.table::fwrite(x = dts_ftd_CPM, file = paste0("03-expression/merged/filtering/UMIcount_filtered_CPM.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_filtered_CPM.txt", " > ", "03-expression/merged/filtering/UMIcount_filtered_CPM.txt.gz"))
# for quick read
write_feather(x = dts_ftd_CPM, sink = paste0("03-expression/merged/filtering/UMIcount_filtered_CPM.feather"))
write.table(x = rownames(dts_ftd_CPM), file = paste0("03-expression/merged/filtering/UMIcount_filtered_CPM.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# 2. cleanFq stat
mt1_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/01-cleandata/merged/cleanFqStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  return(y)
})
names(mt1_LS) <- NULL
mt1_DF <- do.call("rbind", mt1_LS)
dir.create(path = paste0("01-cleandata/merged"), showWarnings = F, recursive = T)
write.table(x = mt1_DF, file = paste0("01-cleandata/merged/cleanFqStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 3. map stat
mt2_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/02-alignment/merged/mapStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  return(y)
})
names(mt2_LS) <- NULL
mt2_DF <- do.call("rbind", mt2_LS)
dir.create(path = paste0("02-alignment/merged"), showWarnings = F, recursive = T)
write.table(x = mt2_DF, file = paste0("02-alignment/merged/mapStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

mt3_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/02-alignment/merged/readsDistri.txt"), header = F, sep = "\t", stringsAsFactors = F)
  return(y)
})
names(mt3_LS) <- NULL
mt3_DF <- do.call("rbind", mt3_LS)
dir.create(path = paste0("02-alignment/merged"), showWarnings = F, recursive = T)
write.table(x = mt3_DF, file = paste0("02-alignment/merged/readsDistri.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 4. expr stat
mt4_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/03-expression/merged/exprStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  return(y)
})
names(mt4_LS) <- NULL
mt4_DF <- do.call("rbind", mt4_LS)
dir.create(path = paste0("03-expression/merged"), showWarnings = F, recursive = T)
write.table(x = mt4_DF, file = paste0("03-expression/merged/exprStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 5. ident
mt5_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/03-expression/merged/cellCluster/Seurat_metaData_pooled.txt"), header = T, sep = "\t", stringsAsFactors = F)
  return(y)
})
names(mt5_LS) <- NULL
mt5_DF <- do.call("rbind", mt5_LS)
dir.create(path = paste0("03-expression/merged/cellCluster"), showWarnings = F, recursive = T)
write.table(x = mt5_DF, file = paste0("03-expression/merged/cellCluster/Seurat_metaData_pooled.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# 6. marker genes
# mt6_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
#   #print(x)
#   y <- read.table(paste0("../", x[1], "/03-expression/merged/cellCluster/Seurat_markerGenes_pooled.txt"), header = T, sep = "\t", stringsAsFactors = F)
#   y$tissue <- gsub("_", " ", gsub("^[0-9][0-9]_", "", x[1]))
#   return(y)
# })
# names(mt6_LS) <- NULL
# mt6_DF <- do.call("rbind", mt6_LS)
# dir.create(path = paste0("03-expression/merged/cellCluster"), showWarnings = F, recursive = T)
# write.table(x = mt6_DF, file = paste0("03-expression/merged/cellCluster/Seurat_markerGenes_pooled.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# 7. cell type order and color
mt7_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../", x[1], "/03-expression/merged/cellCluster/Seurat_cellType_color.txt"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  return(y)
})
names(mt7_LS) <- NULL
mt7_DF <- do.call("rbind", mt7_LS)
dir.create(path = paste0("03-expression/merged/cellCluster"), showWarnings = F, recursive = T)
write.table(x = mt7_DF, file = paste0("03-expression/merged/cellCluster/Seurat_cellType_color.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# 8. cell meta
mt8_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y1 <- read.table(paste0("../", x[1], "/cell_metatable.txt"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  y1$stage <- ifelse(y1$stage == "20w", "19-22w", "11-14w")
  y2 <- read.table(paste0("../", x[1], "/cell_metatable_filtered.txt"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  y2$stage <- ifelse(y2$stage == "20w", "19-22w", "11-14w")
  y3 <- read.table(paste0("../", x[1], "/cell_metatable_filtered_plus.txt"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  y3$stage <- ifelse(y3$stage == "20w", "19-22w", "11-14w")
  y4 <- read.table(paste0("../", x[1], "/cell_metatable_filtered_aligned.txt"), header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  y4$stage <- ifelse(y4$stage == "20w", "19-22w", "11-14w")
  y <- list(y1, y2, y3, y4)
  return(y)
})
names(mt8_LS) <- NULL
# full
cellMeta_final <- do.call("rbind", lapply(mt8_LS, "[[", 1))
write.table(x = cellMeta_final, file = "cell_metatable.txt", row.names = F, col.names = T, quote = F, sep = "\t")
# filtered
cellMeta_final_filtered <- do.call("rbind", lapply(mt8_LS, "[[", 2))
write.table(x = cellMeta_final_filtered, file = "cell_metatable_filtered.txt", row.names = F, col.names = T, quote = F, sep = "\t")
# with cell group
cellMeta_final_filtered_plus <- do.call("rbind", lapply(mt8_LS, "[[", 3))
write.table(x = cellMeta_final_filtered_plus, file = "cell_metatable_filtered_plus.txt", row.names = F, col.names = T, quote = F, sep = "\t")
# with aligned pos
cellMeta_final_filtered_ali <- do.call("rbind", lapply(mt8_LS, "[[", 4))
write.table(x = cellMeta_final_filtered_ali, file = "cell_metatable_filtered_aligned.txt", row.names = F, col.names = T, quote = F, sep = "\t")

# 9. cell type meta
ts_ordered <- read.table("../../pooled_data/All/tissue_ordered.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(ts_ordered)
colnames(ts_ordered) <- "tissue"
ts_ordered$tissue <- tolower(ts_ordered$tissue)
ctMeta_LS <- lapply(mt7_LS, function(x) {
  stages <- unique(x$stage)
  if(length(stages) > 1) {
    y <- subset(x, stage == "14w")
  } else {
    y <- x
  }
  return(y)
})
ctMeta_DF <- do.call("rbind", ctMeta_LS[match(gsub(" ", "_", ts_ordered$tissue), gsub("^[0-9][0-9]_", "", dts$tissue))])[, 1:3]
rownames(ctMeta_DF) <- NULL
write.table(x = ctMeta_DF, file = "cellType_metatable.txt", row.names = F, col.names = T, quote = F, sep = "\t")
