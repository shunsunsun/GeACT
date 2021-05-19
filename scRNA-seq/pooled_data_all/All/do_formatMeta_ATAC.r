
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/All/")

cellMeta_RNA <- read.table(file = "cell_metatable_RNA_global.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_RNA)
cellMeta_ATAC <- read.table(file = "cell_metatable_ATAC.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(cellMeta_ATAC)

# check ident
setdiff(unique(cellMeta_ATAC$ident), unique(cellMeta_RNA$ident))
lapply(unique(cellMeta_ATAC$tissue), function(x) {
  #cat(">", x, "\n")
  y <- setdiff(unique(subset(cellMeta_ATAC, tissue == x, "ident", drop = T)), unique(subset(cellMeta_RNA, tissue == x, "ident", drop = T)))
  return(y)
})
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

# fix embedding
cellMeta_ATAC_fmt_gd <- subset(cellMeta_ATAC_fmt, ! is.na(tSNE_1_ali))
cellMeta_ATAC_fmt_na <- subset(cellMeta_ATAC_fmt, is.na(tSNE_1_ali))
cellMeta_ATAC_fmt_na$tSNE_1_ali <- cellMeta_ATAC_fmt_na$tSNE_1
cellMeta_ATAC_fmt_na$tSNE_2_ali <- cellMeta_ATAC_fmt_na$tSNE_2
cellMeta_ATAC_fmt_na$UMAP_1_ali <- cellMeta_ATAC_fmt_na$UMAP_1
cellMeta_ATAC_fmt_na$UMAP_2_ali <- cellMeta_ATAC_fmt_na$UMAP_2
cellMeta_ATAC_fixed <- rbind(cellMeta_ATAC_fmt_gd, cellMeta_ATAC_fmt_na)
cellMeta_ATAC_fixed <- cellMeta_ATAC_fixed[match(cellMeta_ATAC_fmt$cell, cellMeta_ATAC_fixed$cell), ]
#

data.table::fwrite(x = cellMeta_ATAC_fixed, file = "cell_metatable_ATAC_global.txt", row.names = F, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", "cell_metatable_ATAC_global.txt", " > ", "cell_metatable_ATAC_global.txt.gz"))
