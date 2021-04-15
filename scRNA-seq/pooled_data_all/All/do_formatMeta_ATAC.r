
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

cellMeta_RNA <- read.table(file = "cell_metatable_RNA_global.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_RNA)
cellMeta_ATAC <- read.table(file = "cell_metatable_ATAC.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_ATAC)

# check ident
setdiff(unique(cellMeta_ATAC$ident), unique(cellMeta_RNA$ident))
###
#cellMeta_ATAC[cellMeta_ATAC$group == "TRUE", "group"] <- "T"
###

# change column names
name_DF <- data.frame(old = c("X", "MitoRatio", "BlacklistRatio", "DoubletScore", "nFrags", 
                              "pair_tSNE_1", "pair_tSNE_2", "pair_UMAP_1", "pair_UMAP_2", 
                              "all_tSNE_1", "all_tSNE_2", "all_UMAP_1", "all_UMAP_2"), 
                      new = c("cell", "mitoRatio", "blacklistRatio", "doubletScore", "nFragment", 
                              "tSNE_1_ali", "tSNE_2_ali", "UMAP_1_ali", "UMAP_2_ali", 
                              "tSNE_1_glo", "tSNE_2_glo", "UMAP_1_glo", "UMAP_2_glo"), 
                      stringsAsFactors = F)
colnames(cellMeta_ATAC)[match(name_DF$old, colnames(cellMeta_ATAC))] <- name_DF$new

setdiff(colnames(cellMeta_RNA), colnames(cellMeta_ATAC))
setdiff(colnames(cellMeta_ATAC), colnames(cellMeta_RNA))

col_used <- colnames(cellMeta_RNA)
col_map <- data.frame(ori = setdiff(colnames(cellMeta_RNA), colnames(cellMeta_ATAC)), stringsAsFactors = F)
col_map$new <- c("FRIP", "nPeak", "nFragment", "blacklistRatio", "doubletScore")
col_used[match(col_map$ori, col_used)] <- col_map$new

cellMeta_ATAC_fmt <- cellMeta_ATAC[, col_used]
dim(cellMeta_ATAC_fmt)

data.table::fwrite(x = cellMeta_ATAC_fmt, file = "cell_metatable_ATAC_global.txt", row.names = F, col.names = T, quote = F, sep = "\t", nThread = 10)
