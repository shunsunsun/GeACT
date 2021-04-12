# pool datasets from different stages
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/06_spleen/")

library("parallel")
suppressMessages(library("arrow"))

# options
dts <- data.frame(stgid = c("pooled_data", "pooled_data_14w"), tsid = gsub(".*/", "", getwd()), stringsAsFactors = F)[1, ]
#out <- "X"

# 1.1 UMI all genes (before filtering)
cl <- makeCluster(min(20, nrow(dts)))
dts_LS <- parLapply(cl, split(dts, 1:nrow(dts)), function(x) { 
  print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/03-expression/merged/filtering/UMIcount_unfiltered.txt"), header = T, sep = "\t", 
                  stringsAsFactors = F, row.names = 1, check.names = F, comment.char = "")
  return(y)
})
stopCluster(cl); rm(cl)
names(dts_LS) <- NULL
dts_DF <- do.call("cbind", dts_LS)
dir.create(path = paste0("03-expression/merged/filtering"), showWarnings = F, recursive = T)
data.table::fwrite(x = dts_DF, file = paste0("03-expression/merged/filtering/UMIcount_unfiltered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", "03-expression/merged/filtering/UMIcount_unfiltered.txt", " > ", "03-expression/merged/filtering/UMIcount_unfiltered.txt.gz"))
# for quick read
write_feather(x = dts_DF, sink = paste0("03-expression/merged/filtering/UMIcount_unfiltered.feather"))
write.table(x = rownames(dts_DF), file = paste0("03-expression/merged/filtering/UMIcount_unfiltered.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# 1.2 UMI (after filtering)
# filtering cells
ftc_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/03-expression/merged/filtering/filtering_cells.txt"), header = T, sep = "\t", stringsAsFactors = F)
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
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
# for quick read
write_feather(x = dts_ftd, sink = paste0("03-expression/merged/filtering/UMIcount_filtered.feather"))
write.table(x = rownames(dts_ftd), file = paste0("03-expression/merged/filtering/UMIcount_filtered.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# 2. cleanFq stat
mt1_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/01-cleandata/merged/cleanFqStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
  return(y)
})
names(mt1_LS) <- NULL
mt1_DF <- do.call("rbind", mt1_LS)
dir.create(path = paste0("01-cleandata/merged"), showWarnings = F, recursive = T)
write.table(x = mt1_DF, file = paste0("01-cleandata/merged/cleanFqStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 3. map stat
mt2_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/02-alignment/merged/mapStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
  return(y)
})
names(mt2_LS) <- NULL
mt2_DF <- do.call("rbind", mt2_LS)
dir.create(path = paste0("02-alignment/merged"), showWarnings = F, recursive = T)
write.table(x = mt2_DF, file = paste0("02-alignment/merged/mapStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

mt3_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/02-alignment/merged/readsDistri.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
  return(y)
})
names(mt3_LS) <- NULL
mt3_DF <- do.call("rbind", mt3_LS)
dir.create(path = paste0("02-alignment/merged"), showWarnings = F, recursive = T)
write.table(x = mt3_DF, file = paste0("02-alignment/merged/readsDistri.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 4. expr stat
mt4_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/03-expression/merged/exprStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
  return(y)
})
names(mt4_LS) <- NULL
mt4_DF <- do.call("rbind", mt4_LS)
dir.create(path = paste0("03-expression/merged"), showWarnings = F, recursive = T)
write.table(x = mt4_DF, file = paste0("03-expression/merged/exprStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 5. ident
mt5_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/03-expression/merged/cellCluster/Seurat_metaData.txt"), header = T, sep = "\t", stringsAsFactors = F)
  y <- y[, - grep("^res", colnames(y))]
  u <- read.table(paste0("../../", x[1], "/", x[2], "/03-expression/merged/cellCluster/Seurat_UMAP_embeddings.txt"), header = T, sep = "\t", stringsAsFactors = F)
  colnames(u) <- c("cell", "UMAP_1", "UMAP_2")
  y <- merge(y, u, by = 1, sort = F)
  y <- y[, c(1:13, 23:24, 14:22)]
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
  return(y)
})
names(mt5_LS) <- NULL
mt5_DF <- do.call("rbind", mt5_LS)
dir.create(path = paste0("03-expression/merged/cellCluster"), showWarnings = F, recursive = T)
write.table(x = mt5_DF, file = paste0("03-expression/merged/cellCluster/Seurat_metaData_pooled.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# 6. cell type order and color
mt6_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y <- read.table(paste0("../../", x[1], "/", x[2], "/03-expression/merged/cellCluster/color_DF.txt"), header = F, sep = "\t", stringsAsFactors = F, comment.char = "")[, c(2,1,3)]
  colnames(y) <- c("tissue", "ident", "color")
  y$tissue <- tolower(y$tissue)
  y$stage <- gsub("pooled_data_?", "", as.character(x[1]))
  y$stage[y$stage == ""] <- "20w"
  return(y)
})
names(mt6_LS) <- NULL
mt6_DF <- do.call("rbind", mt6_LS)
dir.create(path = paste0("03-expression/merged/cellCluster"), showWarnings = F, recursive = T)
write.table(x = mt6_DF, file = paste0("03-expression/merged/cellCluster/Seurat_cellType_color.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# 7. cell meta
cellMeta <- merge(ftc_DF[, 1:11], mt1_DF[, c(1,2,5,13)], by.x = "cell", by.y = "V5", sort = F)

sid_info <- read.table("../../datasets/sid_info.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(sid_info)

cellMeta$sid <- gsub("_.*", "", cellMeta$cell)
cellMeta <- merge(cellMeta, sid_info, by = "sid", sort = F)[, -1]
cellMeta <- cellMeta[match(ftc_DF$cell, cellMeta$cell), ]

# add ident, tSNE and UMAP
cellMeta_withIdent <- merge(cellMeta, mt5_DF[, c("cell", "ident", "tSNE_1", "tSNE_2", "UMAP_1", "UMAP_2")], by = "cell", sort = F, all.x = T)
cellMeta_withIdent <- cellMeta_withIdent[match(cellMeta$cell, cellMeta_withIdent$cell), ]
colnames(cellMeta_withIdent)[11:14] <- c("QC", "species", "plate", "seqID")
cellMeta_final <- cellMeta_withIdent[, c(1,14,12,15,16,13,11,17,2:8,10,18,19,20,21)]

# add stage
cellMeta_final <- merge(cellMeta_final, ftc_DF[, c("cell", "stage")], by = "cell", sort = F)

all(colnames(dts_DF) == cellMeta_final$cell)

# full
write.table(x = cellMeta_final, file = "cell_metatable.txt", row.names = F, col.names = T, quote = F, sep = "\t")
# filtered
cellMeta_final_filtered <- subset(cellMeta_final, QC)
write.table(x = cellMeta_final_filtered, file = "cell_metatable_filtered.txt", row.names = F, col.names = T, quote = F, sep = "\t")

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

cellMeta_final_filtered_plus <- cellMeta_final_filtered
cellMeta_final_filtered_plus$group <- ident2clgrp(cellMeta_final_filtered_plus$ident)
cellMeta_final_filtered_plus <- cellMeta_final_filtered_plus[, c(1:8,22,9:21)]
write.table(x = cellMeta_final_filtered_plus, file = "cell_metatable_filtered_plus.txt", row.names = F, col.names = T, quote = F, sep = "\t")

# add aligned pos
ts_id <- Hmisc::capitalize(gsub(".*[0-9][0-9]_", "", getwd()))
if(file.exists(paste0("../../pooled_data_14w/cblast_result/", ts_id, ".csv"))) {
  cat("> Use aligned embeddings.\n")
  cblast_res <- read.table(file = paste0("../../pooled_data_14w/cblast_result/", ts_id, ".csv"), header = T, sep = ",", stringsAsFactors = F, row.names = 1)
  rownames(cblast_res) <- gsub(".", "-", rownames(cblast_res), fixed = T)
  cblast_res_cub <- cblast_res[, c("tSNE1", "tSNE2", "UMAP1", "UMAP2")]
} else {
  cat("> Use original embeddings.\n")
  cblast_res_cub <- cellMeta_final_filtered_plus[, c("tSNE_1", "tSNE_2", "UMAP_1", "UMAP_2")]
  rownames(cblast_res_cub) <- cellMeta_final_filtered_plus$cell
}
colnames(cblast_res_cub) <- c("tSNE_1_ali", "tSNE_2_ali", "UMAP_1_ali", "UMAP_2_ali")
cellMeta_final_filtered_ali <- merge(cellMeta_final_filtered_plus, cblast_res_cub, by.x = "cell", by.y = 0, sort = F)
write.table(x = cellMeta_final_filtered_ali, file = "cell_metatable_filtered_aligned.txt", row.names = F, col.names = T, quote = F, sep = "\t")

# 8. cell type meta
ts_ordered <- read.table("../../pooled_data/All/tissue_ordered.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(ts_ordered)
colnames(ts_ordered) <- "tissue"
ts_ordered$tissue <- tolower(ts_ordered$tissue)
ctMeta <- do.call("rbind", split(mt6_DF, mt6_DF$tissue)[ts_ordered$tissue])
rownames(ctMeta) <- NULL
write.table(x = ctMeta, file = "cellType_metatable.txt", row.names = F, col.names = T, quote = F, sep = "\t")
