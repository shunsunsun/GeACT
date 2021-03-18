# cell meta

setwd("~/lustre/06-Human_cell_atlas/pooled_data/KD/")

ftc_DF <- read.table(file = "03-expression/merged/filtering/filtering_cells.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(ftc_DF)[1] <- "cell"
mt1_DF <- read.table(file = "01-cleandata/merged/cleanFqStat.txt", header = F, sep = "\t", stringsAsFactors = F)
mt5_DF <- read.table(file = "03-expression/merged/cellCluster/Seurat_metaData.txt", header = T, sep = "\t", stringsAsFactors = F)
mt5_DF <- mt5_DF[, ! colnames(mt5_DF) %in% c("res.0.6", "cond", "cond_cluster", "cluster_pure")]
u <- read.table("03-expression/merged/cellCluster/Seurat_UMAP_embeddings.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(u) <- c("cell", "UMAP_1", "UMAP_2")
mt5_DF <- merge(mt5_DF, u, by = 1, sort = F)
mt5_DF <- mt5_DF[, c(1:13, 23:24, 14:22)]

cellMeta <- merge(ftc_DF, mt1_DF[, c(1,2,5,13)], by.x = "cell", by.y = "V5", sort = F)

sid_info <- read.table("../../datasets/sid_info.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(sid_info)
if(! "XC" %in% sid_info$sid) { sid_info <- rbind(sid_info, c("XC", "small intestine", "small intestine"))}

cellMeta$sid <- gsub("_.*", "", cellMeta$cell)
cellMeta <- merge(cellMeta, sid_info, by = "sid", sort = F)[, -1]
cellMeta <- cellMeta[match(ftc_DF$cell, cellMeta$cell), ]

# add ident, tSNE and UMAP
cellMeta_withIdent <- merge(cellMeta, mt5_DF[, c("cell", "ident", "tSNE_1", "tSNE_2", "UMAP_1", "UMAP_2")], by = "cell", sort = F, all.x = T)
cellMeta_withIdent <- cellMeta_withIdent[match(cellMeta$cell, cellMeta_withIdent$cell), ]
colnames(cellMeta_withIdent)[11:14] <- c("QC", "species", "plate", "seqID")
cellMeta_final <- cellMeta_withIdent[, c(1,14,12,15,16,13,11,17,2:8,10,18,19,20,21)]

###
cellMeta_final <- subset(cellMeta_final, ! grepl("NR2F1", cell))
###

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
cellMeta_final_filtered_plus <- cellMeta_final_filtered_plus[, c(1:8,21,9:20)]
write.table(x = cellMeta_final_filtered_plus, file = "cell_metatable_filtered_plus.txt", row.names = F, col.names = T, quote = F, sep = "\t")
