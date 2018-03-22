#!/usr/bin/Rscript

#setwd("/lustre/user/tianf/06-Human_cell_atlas/HCA-web/05_scripts/cell_blast")

# get options
cmds <- commandArgs(trailingOnly = T)
qef <- cmds[1]
qft <- cmds[2]
ref <- cmds[3]
ctf <- cmds[4]
gld <- cmds[5]
mtd <- cmds[6]
odr <- cmds[7]

### just for development
# {
#   qef <- "demo/query.txt"
#   qft <- "ensembl"
#   ref <- "demo/reference.txt"
#   ctf <- "null"
#   gld <- "demo/signature.txt"
#   mtd <- "PCC"
#   odr <- "test"
# }
###

do_norm <- function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
}

do_highlyVar <- function(m) {
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) { data.frame(bin_median=median(x$dispersion), bin_mad=mad(x$dispersion)) })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}

query <- read.table(file = qef, header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)  # genes * cells
reference <- read.table(file = ref, header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)  # genes * cell types
dim(query)
dim(reference)

# parse query
query_t <- t(query) # cells * genes
query_normed <- do_norm(query_t)
# extract signature genes
if(gld != "null") {
  cat("[info] use the signature genes from users.\n")
  signature_genes <- read.table(file = gld, header = F, sep = "\t", stringsAsFactors = F)[,1]
  signature_genes <- unique(signature_genes)
} else {
  cat("[info] use the signature genes from 1000 highest variable genes in query.\n")
  query_var_DF <- do_highlyVar(query_normed)
  high_var_genes <- rownames(query_var_DF)[order(query_var_DF$dispersion_norm, decreasing = T)][1:1000]
  signature_genes <- high_var_genes
}
signature_genes_valid <- signature_genes[signature_genes%in%colnames(query_normed) & signature_genes%in%rownames(reference)]
query_normed_signature <- t(query_normed[, signature_genes_valid]) # genes * cells
reference_signature <- reference[signature_genes_valid, ] # genes * cells

# calculate similarity
if(mtd=="PCC") {
  cat("[info] use the maxmum Spearman correlation to identify cell types.\n")
  cor_MT <- cor(query_normed_signature, reference_signature, method = "spearman")
  cellToType <- t(apply(cor_MT, 1, function(x) {
    ya <- colnames(cor_MT)[which.max(x)]
    yb <- max(x)
    yc <- sort(x, decreasing = T)[2]
    yd <- mean(sort(x, decreasing = T)[-1])
    return(c(ya, yb, yc, yd))
  }))
  cellToType <- data.frame(cell=rownames(cellToType), type=cellToType[,1],
                           similarity=as.numeric(cellToType[,2]), nextValue=as.numeric(cellToType[,3]), negative=as.numeric(cellToType[,4]), stringsAsFactors = F)
}

cellToType_tmp <- cellToType
cellToType_tmp[,3:5] <- format(cellToType_tmp[,3:5], digits = 4, scientific = F)
write.table(x = cellToType_tmp, file = paste0(odr,"/cellToType.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
