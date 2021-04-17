# average expression
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

suppressMessages(library("arrow"))

samplingPos <- "."
OUT <- paste0("03-expression/merged/filtering/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

# 1. using the genes before filtering ----
# Load the expr matrix (CPM)
expr_data_normed <- read_feather(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_cellFiltered_CPM.feather"))
expr_data_normed <- as.data.frame(expr_data_normed)
expr_data_normed_gene <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_cellFiltered_CPM.gene"), header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_data_normed) <- expr_data_normed_gene$V1

# Load the cell metatable
cellMetaData <- read.table(file = "../../pooled_data_all/All/cell_metatable_RNA_global.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
cellMetaData <- cellMetaData[colnames(expr_data_normed), ]
all(colnames(expr_data_normed) == rownames(cellMetaData))
cellMetaData$ts_ident <- paste(cellMetaData$tissue, cellMetaData$ident, sep = ".")
#cellMetaData$ts_clgrp <- Hmisc::capitalize(paste(cellMetaData$tissue, ident2clgrp(cellMetaData$ident), sep = "."))
#length(unique(cellMetaData$ts_clgrp))
#View(unique(cellMetaData$ts_clgrp))

# avg
expr_data_avg <- sapply(split(rownames(cellMetaData), cellMetaData$ts_ident), function(x) { y <- rowMeans(expr_data_normed[, x, drop = F]) })
expr_data_avg <- as.data.frame(expr_data_avg)

data.table::fwrite(x = expr_data_avg, file = paste0(OUT, "/UMIcount_cellFiltered_avgCPM_byCt.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", OUT, "/UMIcount_cellFiltered_avgCPM_byCt.txt", " > ", OUT, "/UMIcount_cellFiltered_avgCPM_byCt.txt.gz"))
# for quick read
write_feather(x = expr_data_avg, sink = paste0(OUT, "/UMIcount_cellFiltered_avgCPM_byCt.feather"))
write.table(x = rownames(expr_data_avg), file = paste0(OUT, "/UMIcount_cellFiltered_avgCPM_byCt.gene"), row.names = F, col.names = F, quote = F, sep = "\t")

# 2. using the genes after filtering ----
# Load the expr matrix (CPM)
expr_data_normed <- read_feather(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered_CPM.feather"))
expr_data_normed <- as.data.frame(expr_data_normed)
expr_data_normed_gene <- read.table(paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered_CPM.gene"), header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_data_normed) <- expr_data_normed_gene$V1

# Load the cell metatable
cellMetaData <- read.table(file = "../../pooled_data_all/All/cell_metatable_RNA_global.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
cellMetaData <- cellMetaData[colnames(expr_data_normed), ]
all(colnames(expr_data_normed) == rownames(cellMetaData))
cellMetaData$ts_ident <- paste(cellMetaData$tissue, cellMetaData$ident, sep = ".")
#cellMetaData$ts_clgrp <- Hmisc::capitalize(paste(cellMetaData$tissue, ident2clgrp(cellMetaData$ident), sep = "."))
#length(unique(cellMetaData$ts_clgrp))
#View(unique(cellMetaData$ts_clgrp))

# avg
expr_data_avg <- sapply(split(rownames(cellMetaData), cellMetaData$ts_ident), function(x) { y <- rowMeans(expr_data_normed[, x, drop = F]) })
expr_data_avg <- as.data.frame(expr_data_avg)

data.table::fwrite(x = expr_data_avg, file = paste0(OUT, "/UMIcount_filtered_avgCPM_byCt.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", OUT, "/UMIcount_filtered_avgCPM_byCt.txt", " > ", OUT, "/UMIcount_filtered_avgCPM_byCt.txt.gz"))
# for quick read
write_feather(x = expr_data_avg, sink = paste0(OUT, "/UMIcount_filtered_avgCPM_byCt.feather"))
write.table(x = rownames(expr_data_avg), file = paste0(OUT, "/UMIcount_filtered_avgCPM_byCt.gene"), row.names = F, col.names = F, quote = F, sep = "\t")
