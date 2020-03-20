# pool datasets from different seqID
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("parallel")

# options
dts <- data.frame(tissue = list.files(path = "..", pattern = "^[0-9]"), stringsAsFactors = F)
#out <- "X"

# 1.1 UMI (before filtering)
cl <- makeCluster(min(nrow(dts), 10))
dts_LS <- parLapply(cl, split(dts, 1:nrow(dts)), function(x) { 
  print(x)
  y <- read.table(paste0("../", x[1], "/03-expression/merged/filtering/UMIcount_unfiltered.txt"), header = T, sep = "\t", 
                  stringsAsFactors = F, row.names = 1, check.names = F, comment.char = "")
  return(y)
})
stopCluster(cl); rm(cl)
names(dts_LS) <- NULL
dts_DF <- do.call("cbind", dts_LS)
dir.create(path = paste0("03-expression/merged/filtering"), showWarnings = F, recursive = T)
# write.table vs fwrite (1000 rows * 31392 columns: 186.282s ~ 4.647s)
data.table::fwrite(x = dts_DF, file = paste0("03-expression/merged/filtering/UMIcount_unfiltered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_unfiltered.txt", " > ", "03-expression/merged/filtering/UMIcount_unfiltered.txt.gz"))

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

# filtering genes
nCell_expressed <- rowSums(dts_cellftd > 0)
dts_ftd <- dts_cellftd[nCell_expressed >= 10, ]
dim(dts_ftd)
data.table::fwrite(x = dts_ftd, file = paste0("03-expression/merged/filtering/UMIcount_filtered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_filtered.txt", " > ", "03-expression/merged/filtering/UMIcount_filtered.txt.gz"))

# norm
dts_ftd_CPM <- sweep(x = dts_ftd, MARGIN = 2, STATS = colSums(dts_ftd), FUN = "/") * 1e6
data.table::fwrite(x = dts_ftd_CPM, file = paste0("03-expression/merged/filtering/UMIcount_filtered_CPM.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 15)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_filtered_CPM.txt", " > ", "03-expression/merged/filtering/UMIcount_filtered_CPM.txt.gz"))

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
  y <- read.table(paste0("../", x[1], "/03-expression/merged/cellCluster/Seurat_metaData.txt"), header = T, sep = "\t", stringsAsFactors = F)
  y <- y[, - grep("^res", colnames(y))]
  return(y)
})
names(mt5_LS) <- NULL
mt5_DF <- do.call("rbind", mt5_LS)
dir.create(path = paste0("03-expression/merged/cellCluster"), showWarnings = F, recursive = T)
write.table(x = mt5_DF, file = paste0("03-expression/merged/cellCluster/Seurat_metaData_pooled.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# 6. meta
cellMeta <- merge(ftc_DF, mt1_DF[, c(1,2,5,13)], by.x = "cell", by.y = "V5", sort = F)

### meta from Github
# meta_table <- read.table(file = "meta_table.txt", header = T, sep = "\t", stringsAsFactors = F)
# dim(meta_table)
# # only use HCA
# meta_table <- meta_table[grep("_HCA_", meta_table$plate), ]
# meta_table <- meta_table[, c("plate", "tissue", "samplingPos")]
# meta_table <- meta_table[! duplicated(meta_table), ]
# # rename
# meta_table$plate <- gsub("-CGM-", "_CGM_", meta_table$plate)
# meta_table$plate <- gsub("_CGM_1-", "_C-A", meta_table$plate)
# meta_table$plate <- gsub("_CGM_2-", "_C-B", meta_table$plate)
# meta_table$plate <- gsub("_CGM_3-", "_C-C", meta_table$plate)
# meta_table$plate <- gsub("_CGM_4-", "_C-D", meta_table$plate)
# 
# setdiff(meta_table$plate, cellMeta$V2)
# setdiff(cellMeta$V2, meta_table$plate)
# cellMeta <- merge(cellMeta, meta_table, by.x = "V2", by.y = "plate", sort = F)
###

sid_info <- read.table("../../datasets/sid_info.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(sid_info)

cellMeta$sid <- gsub("_.*", "", cellMeta$cell)
cellMeta <- merge(cellMeta, sid_info, by = "sid", sort = F)[, -1]
cellMeta <- cellMeta[match(ftc_DF$cell, cellMeta$cell), ]

# add ident and tSNE
cellMeta_withIdent <- merge(cellMeta, mt5_DF[, c("cell", "ident", "tSNE_1", "tSNE_2")], by = "cell", sort = F, all.x = T)
cellMeta_withIdent <- cellMeta_withIdent[match(cellMeta$cell, cellMeta_withIdent$cell), ]
colnames(cellMeta_withIdent)[11:14] <- c("QC", "species", "plate", "seqID")
cellMeta_final <- cellMeta_withIdent[, c(1,14,12,15,16,13,11,17,2:8,10,18,19)]

all(colnames(dts_DF) == cellMeta_final$cell)

# full
write.table(x = cellMeta_final, file = "cell_metatable.txt", row.names = F, col.names = T, quote = F, sep = "\t")
# filtered
cellMeta_final_filtered <- subset(cellMeta_final, QC)
write.table(x = cellMeta_final_filtered, file = "cell_metatable_filtered.txt", row.names = F, col.names = T, quote = F, sep = "\t")
