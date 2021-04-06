# harmony
suppressMessages(library("arrow"))
suppressMessages(library("harmony"))
library("Matrix")

# load expression matrix
expr_ftd <- read_feather("03-expression/merged/filtering/UMIcount_filtered.feather")
expr_ftd <- as.data.frame(expr_ftd)
expr_ftd_gene <- read.table("03-expression/merged/filtering/UMIcount_filtered.gene", header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_ftd) <- expr_ftd_gene$V1

# load cell meta
cellMeta <- read.table("cell_metatable_filtered_plus.txt", header = T, sep = "\t", stringsAsFactors = F)

# norm
expr_norm <- sweep(expr_ftd, 2, colSums(expr_ftd), "/") * 1e6
expr_norm_MT <- as.matrix(expr_norm)

# alignment
harmony_embeddings <- HarmonyMatrix(expr_norm_MT, cellMeta, "stage")
harmony_embeddings_cub <- as.data.frame(harmony_embeddings[, 1:2])
