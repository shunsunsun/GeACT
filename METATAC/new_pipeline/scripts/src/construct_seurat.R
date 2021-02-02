#!/usr/bin/env r
suppressMessages({
  library(docopt)
})


# Argument parsing -------------------------------------------------------------
doc <- "Usage: construct_seurat [-h] [-s STAGE] [-t TISSUE] [--root ROOT] [--buildstage BUILDSTAGE] [--maxcutoff MAXCO]

-h --help           show help message
-s --stage STAGE    stage of sample [default: 19-22w]
-t --tissue TISSUE  tissue name [default: 02_small_intestine]
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
--buildstage BUILDSTAGE target stage to build [default: 19-22w]
--maxcutoff MAXCO   Max cutoff when finding variable features [default: 8]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

suppressMessages({
  library(Seurat)
  library(tidyverse)
})

# envs -------------------------------------------------------------------------
stage <- opt$stage
tissue <- opt$tissue
root <- opt$root
buildstage <- opt$buildstage
max_cutoff <- as.integer(opt$maxcutoff)

.data <- "data"
utils_file <- paste(root, "scripts/src/utils.R", sep = "/")

organ_wd <- paste(root, .data, stage, tissue, sep = "/")
rna_wd <- paste(organ_wd, "RNA", buildstage, sep = "/")

source(utils_file)
setwd(rna_wd)

filtered_cell_file <- "UMIcount_cellFiltered.txt"
cell_meta_file <- "Seurat_metaData.txt"

expr_filtered_cell <- read.csv(filtered_cell_file, header = T, sep = "\t", row.names = 1, check.names = F)

cell_meta <- read_tsv(cell_meta_file, col_names = T)
cell_meta <- column_to_rownames(cell_meta, "cell")[colnames(expr_filtered_cell), ]

ts <- CreateSeuratObject(expr_filtered_cell, project = paste(stage, tissue, sep = "_"), names.field = 0, meta.data = cell_meta, min.cells = 10, min.features = 0)
Idents(ts) <- ts[["ident"]]

ts <- NormalizeData(ts)
pdf("VF.pdf", width = 5, height = 5)
print(ggplot(data = rowMeans(ts, slot = "counts") %>% my_clip(lb = 0.1, ub = max_cutoff) %>% as_tibble, aes(x=value)) + geom_histogram(bins = 20))
dev.off()
ts <- FindVariableFeatures(ts, selection.method = "vst", mean.cutoff = c(0.1, max_cutoff))

saveRDS(ts, "expr_RNA.rds")
