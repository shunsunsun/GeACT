
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

cellMeta_RNA <- read.table(file = "cell_metatable_filtered_aligned.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_RNA)
cellMeta_ATAC <- read.table(file = "cell_metatable_ATAC.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_ATAC)
###
cellMeta_ATAC[cellMeta_ATAC$group == "TRUE", "group"] <- "T"
###
colnames(cellMeta_ATAC)[1] <- "cell"
colnames(cellMeta_ATAC)[colnames(cellMeta_ATAC) == "MitoRatio"] <- "mitoRatio"
colnames(cellMeta_ATAC)[colnames(cellMeta_ATAC) == "BlacklistRatio"] <- "blacklistRatio"
colnames(cellMeta_ATAC)[colnames(cellMeta_ATAC) == "DoubletScore"] <- "doubletScore"
colnames(cellMeta_ATAC)[colnames(cellMeta_ATAC) == "nFrags"] <- "nFragment"
cellMeta_ATAC$tSNE_1_ali <- cellMeta_ATAC$tSNE_1
cellMeta_ATAC$tSNE_2_ali <- cellMeta_ATAC$tSNE_2
cellMeta_ATAC$UMAP_1_ali <- cellMeta_ATAC$UMAP_1
cellMeta_ATAC$UMAP_2_ali <- cellMeta_ATAC$UMAP_2

setdiff(colnames(cellMeta_RNA), colnames(cellMeta_ATAC))
setdiff(colnames(cellMeta_ATAC), colnames(cellMeta_RNA))

col_used <- colnames(cellMeta_RNA)
col_map <- data.frame(ori = setdiff(colnames(cellMeta_RNA), colnames(cellMeta_ATAC)), stringsAsFactors = F)
col_map$new <- c("FRIP", "nPeak", "nFragment", "blacklistRatio", "doubletScore")
col_used[match(col_map$ori, col_used)] <- col_map$new
col_used <- c(col_used, c("all_tSNE_1", "all_tSNE_2", "all_UMAP_1", "all_UMAP_2"))

cellMeta_ATAC_new <- cellMeta_ATAC[, col_used]
dim(cellMeta_ATAC_new)
colnames(cellMeta_ATAC_new)[27:30] <- c("tSNE_1_glo", "tSNE_2_glo", "UMAP_1_glo", "UMAP_2_glo")

data.table::fwrite(x = cellMeta_ATAC_new, file = "cell_metatable_ATAC_global.txt", row.names = F, col.names = T, quote = F, sep = "\t", nThread = 10)
