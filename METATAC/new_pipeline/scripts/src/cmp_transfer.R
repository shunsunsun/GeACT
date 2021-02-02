#!/usr/bin/env r
suppressMessages({
  library(docopt)
})


# Argument parsing -------------------------------------------------------------
doc <- "Usage: cmp_tranfer.R [-h] [-t TISSUE] [--root ROOT]

-h --help           show help message
-t --tissue TISSUE  tissue name [default: 02_small_intestine]
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

suppressMessages({
  library(tidyverse)
})

stage <- "11-14w"
tissue <- opt$tissue
root <- opt$root

.data <- "data"

organ_wd <- paste(root, .data, stage, tissue, sep = "/")
results_wd <- paste(organ_wd, "results", sep = "/")

utils_file <- paste(root, "scripts/src/utils.R", sep = "/")
source(utils_file)

setwd(results_wd)

tr_to_11_14w <- read.table("filtered_cellMeta.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "") %>% column_to_rownames(var = "X")
tr_to_19_22w <- read.table("filtered_cellMeta_A.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "") %>% column_to_rownames(var = "X")

dir.create("cmp_transfer_stage", showWarnings = F)

setwd("cmp_transfer_stage")

lev <- unique(c(tr_to_11_14w$ident, tr_to_19_22w$ident))
p <- cmp_transfer(tr_to_11_14w$ident, tr_to_19_22w$ident, lev = lev, name1 = "11_14w", name2 = "19-22w")

saveNetwork(p, "cmp_ident.html")

lev_g <- unique(c(tr_to_11_14w$group, tr_to_19_22w$group))
p_g <- cmp_transfer(tr_to_11_14w$group, tr_to_19_22w$group, lev = lev_g, name1 = "11_14w", name2 = "19-22w")

saveNetwork(p_g, "cmp_group.html")
