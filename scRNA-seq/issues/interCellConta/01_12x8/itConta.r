
#setwd("~/lustre/06-Human_cell_atlas/datasets/F1_DT1_45/")
OUT <- "03-expression/itConta"

library("data.table")
suppressMessages(library("reshape2"))
library("parallel")

batch_id <- commandArgs(trailingOnly = T)
if(length(batch_id) != 1) {
  cat("Remove UMI contamination for cells in the same plate.\n")
  cat("Usage:   Rscript do_itConta.r plate_id\n")
  cat("Example: Rscript do_itConta.r SD-1\n")
  q("no")
}

dir.create(path = paste0(OUT, "/", batch_id), showWarnings = F, recursive = T)

tabFiles <- system("ls 03-expression/human/*/*.tab | sort -V", intern = TRUE)
batches <- gsub(".*/", "", gsub("_cell.*", "", tabFiles))
cells <- gsub(".*/", "", gsub("_htseq.*", "", tabFiles))
tab_DF <- data.frame(batch = batches, cell = cells, fl = tabFiles, stringsAsFactors = F)
tab_DF_sub <- subset(tab_DF, batch == batch_id)

### 1. tab to uc ----
cl <- makeCluster(20)
clusterExport(cl, c("fread", "acast"))

input_fl <- tab_DF_sub$fl
cat("Cell number in this plate:", length(input_fl), "\n")

cat("[1.1] Merge tab files\n")
tb_LS <- parLapply(cl, input_fl, function(x) {  # for each cell
  cat(">", x, "\n")
  ka <- fread(file = x, header = F, sep = "\t", stringsAsFactors = F, verbose = F)
  ka <- ka[grep("^__", ka$V1, invert = T), ]
  if(nrow(ka) > 0) {
    ka <- data.frame(cell = gsub(".*/", "", gsub("_htseq.tab", "", x)), ka, stringsAsFactors = F)
  } else {
    cat("No lines in the tab file.\n")
  }
})
tb_DF <- do.call("rbind", tb_LS)
colnames(tb_DF) <- c("cell", "gene", "UMI", "num")

cat("[1.2] Calculate umi_by_cell for each gene\n")
tb_DF_LS <- split(tb_DF[, c("UMI", "cell", "num")], tb_DF$gene)
umiCell_LS <- parLapply(cl, tb_DF_LS, function(xi) {
  y <- acast(data = xi, formula = UMI ~ cell, value.var = "num", fill = 0)
})
stopCluster(cl); rm(cl)

saveRDS(object = umiCell_LS, file = paste0(OUT, "/", batch_id, "/", "umiCell_LS.rds"))
#umiCell_LS <- readRDS(file = paste0(OUT, "/", batch_id, "/", "umiCell_LS.rds"))
rm(tb_LS, tb_DF, tb_DF_LS)
invisible(gc())

### 2. de conta ----
cutoff <- 0.75
rescue <- F

cl <- makeCluster(20)
clusterExport(cl, c("cutoff", "rescue", "acast"))

cat("[2] mark contamination\n")
deConta_res <- parLapply(cl, umiCell_LS, function(x) {  # for each gene
  ### x <- umiCell_LS[[1000]]
  x1 <- x[apply(x>0, 1, sum) > 1, , drop = F]
  if(nrow(x1) == 0) {
    #cat("No UMI exist in more than 1 cells, skip.\n")
    x_res <- x
    stat_res_fs <- data.frame(row.names = colnames(x), 
                              ml_contaRow = rep(0, ncol(x)), ml_contaCol = rep(0, ncol(x)), ml_contaOth = rep(0, ncol(x)), 
                              tp_contaRow = rep(0, ncol(x)), tp_contaCol = rep(0, ncol(x)), tp_contaOth = rep(0, ncol(x)))
  } else {
    x1_res_LS <- apply(x1, 1, function(xi) {  # for each UMI
      ###xi <- x1[1, ]
      ct_max <- max(xi)
      ct_max_ele <- xi[xi == ct_max]
      #print(ct_max_ele)
      xs <- xi
      if(length(ct_max_ele) > 1) {
        #cat("More than 1 cells show max UMI count, skip.\n")
      } else {
        ct_ratio <- xi / ct_max
        xs[ct_ratio > 0 & ct_ratio < cutoff] <- 0
        if(rescue) {
          ct_max_ind <- which.max(xi)
          xs[ct_max_ind] <- xs[ct_max_ind] + sum(xi - xs)
        }
      }
      # conta stat
      diff_DF <- data.frame(cell = names(xi), xi, xs, diff = xi - xs, row.names = NULL, stringsAsFactors = F)
      diff_DF <- subset(diff_DF, diff > 0)
      if(nrow(diff_DF) > 0) {
        diff_DF$outer_id <- ceiling(as.numeric(gsub(".*cell", "", diff_DF$cell)) / 8)
        diff_DF$inner_id <- as.numeric(gsub(".*cell", "", diff_DF$cell)) %% 8; diff_DF$inner_id[diff_DF$inner_id == 0] <- 8
        diff_DF$outer_mx <- ceiling(as.numeric(gsub(".*cell", "", names(ct_max_ele))) / 8)
        diff_DF$inner_mx <- as.numeric(gsub(".*cell", "", names(ct_max_ele))) %% 8; diff_DF$inner_mx[diff_DF$inner_mx == 0] <- 8
        diff_DF$conta <- apply(diff_DF[, -1], 1, function(x) {
          if(x[4]!=x[6] && x[5]!=x[7]) { y = "contaOth" }
          else if(x[5]==x[7]) { y = "contaRow" }  # same inner
          else if(x[4]==x[6]) { y = "contaCol" }  # same outer
          return(y)
        })
      }
      return(list(xs, diff_DF))
    })  # after loop for UMI
    
    x1_res <- t(sapply(x1_res_LS, "[[", 1)) # combine all UMI (result: umi by cell)
    x_res <- rbind(x1_res, x[apply(x>0, 1, sum) == 1, ])
    
    stat_res <- do.call("rbind", lapply(x1_res_LS, "[[", 2))  # combine all UMI
    stat_res$UMI <- gsub("\\.[0-9]+", "", rownames(stat_res))
    if(nrow(stat_res) > 0) {
      stat_res$conta <- factor(stat_res$conta, levels = c("contaRow", "contaCol", "contaOth"))
      # molecule level
      stat_res_fm <- acast(stat_res, formula = cell ~ conta, value.var = "diff", fun.aggregate = sum, drop = F)
      # type level
      stat_res_ft <- acast(stat_res, formula = cell ~ conta, value.var = "UMI", fun.aggregate = function(x) { length(unique(x)) }, drop = F)
    } else {
      stat_res_fm <- data.frame(row.names = colnames(x), 
                                contaRow = rep(0, ncol(x)), contaCol = rep(0, ncol(x)), contaOth = rep(0, ncol(x)), stringsAsFactors = F)
      stat_res_ft <- data.frame(row.names = colnames(x), 
                                contaRow = rep(0, ncol(x)), contaCol = rep(0, ncol(x)), contaOth = rep(0, ncol(x)), stringsAsFactors = F)
    }
    stat_res_fs <- cbind(stat_res_fm, stat_res_ft)
    colnames(stat_res_fs) <- c("ml_contaRow", "ml_contaCol", "ml_contaOth", "tp_contaRow", "tp_contaCol", "tp_contaOth")
  }
  # add origin
  stat_res_fi <- merge(data.frame(ml_total = colSums(x), tp_total = colSums(x>0)), stat_res_fs, by = 0, sort = F, all.x = T)
  stat_res_fi[is.na(stat_res_fi)] <- 0
  stat_res_fi$ml_origin <- stat_res_fi$ml_total - rowSums(stat_res_fi[, c("ml_contaRow", "ml_contaCol", "ml_contaOth")])
  stat_res_fi$tp_origin <- stat_res_fi$tp_total - rowSums(stat_res_fi[, c("tp_contaRow", "tp_contaCol", "tp_contaOth")])
  colnames(stat_res_fi)[1] <- "cell"
  stat_res_fi <- stat_res_fi[, c("cell", 
                                 "ml_total", "ml_contaRow", "ml_contaCol", "ml_contaOth", "ml_origin", 
                                 "tp_total", "tp_contaRow", "tp_contaCol", "tp_contaOth", "tp_origin")]
  return(list(x_res, stat_res_fi))
})  # after loop for gene
stopCluster(cl); rm(cl)

saveRDS(object = deConta_res, file = paste0(OUT, "/", batch_id, "/", "deConta_res.rds"))
#deConta_res <- readRDS(file = paste0(OUT, "/", batch_id, "/", "deConta_res.rds"))
rm(umiCell_LS, cutoff, rescue)
invisible(gc())

### 3. conta_res to tab ----
cl <- makeCluster(20)

cat("[3.1] calc expression matrix\n")
deConta_expr_LS <- lapply(tab_DF_sub$cell, function(cs) { # for each cell
  ### cs <- tab_DF_sub$cell[1]
  cat(">", cs, "\n")
  ys <- parLapply(cl, deConta_res, function(x) {  # for each gene
    ### x <- deConta_res[[1]]
    xx <- x[[1]]
    if(cs %in% colnames(xx)) {
      y <- xx[, cs, drop = F]
      y <- y[y > 0, , drop = F]
      y <- data.frame(UMI = rownames(y), num = y[, 1], stringsAsFactors = F, check.names = F)
    } else {
      y <- data.frame(UMI = character(), num = integer(), stringsAsFactors = F, check.names = F)
    }
    return(y)
  })
  cat("Merge expression for all genes...\n")
  ys <- ys[sapply(ys, nrow) > 0]
  if(length(ys) > 0) {
    yt <- melt(ys, measure.vars = NULL)
    colnames(yt) <- c("UMI", "num", "gene")
    yt <- yt[, c("gene", "UMI", "num")]
  } else {
    yt <- data.frame(gene = character(), UMI = character(), num = integer(), stringsAsFactors = F)
  }
  return(yt)
})
stopCluster(cl); rm(cl)
names(deConta_expr_LS) <- tab_DF_sub$cell

saveRDS(object = deConta_expr_LS, file = paste0(OUT, "/", batch_id, "/", "deConta_expr_LS.rds"))
#deConta_expr_LS <- readRDS(file = paste0(OUT, "/", batch_id, "/", "deConta_expr_LS.rds"))

cat("[3.2] write files\n")
dir.create(path = paste0(OUT, "/", batch_id, "/exprAfDC"), showWarnings = F, recursive = T)
for(i in seq_along(deConta_expr_LS)) {
  cs <- names(deConta_expr_LS)[i]
  #cat(">", cs, "\n")
  write.table(x = deConta_expr_LS[[i]], file = paste0(OUT, "/", batch_id, "/exprAfDC/", cs, "_htseq.tab"), 
              row.names = F, col.names = F, quote = F, sep = "\t")
}

### 4. UMI collapse ----
geneID <- read.table("../../Genomes/human/genes.txt", header = F, sep = "\t", stringsAsFactors = F)

cl <- makeCluster(20)
clusterExport(cl, c("OUT", "batch_id", "geneID"))

cat("[4] UMI collapsing...\n")
dir.create(path = paste0(OUT, "/", batch_id, "/exprRmdup"), showWarnings = F, recursive = T)
expr_rmdup_ensmebl <- parSapply(cl, names(deConta_expr_LS), function(cs) {
  ### cs <- names(deConta_expr_LS)[1]
  # 1
  outfile <- paste0(OUT, "/", batch_id, "/exprRmdup/", cs, "_UMIcount.txt")
  system(paste0("perl scripts/UMIcollapse.pl ", OUT, "/", batch_id, "/exprAfDC/", cs, "_htseq.tab ", 4, " > ", outfile))
  # 2
  umicount <- tryCatch(read.table(outfile, header = F, sep = "\t", stringsAsFactors = F), error=function(e) data.frame(x = integer()))
  if(nrow(umicount) > 0) {
    expr_allGenes <- merge(geneID, umicount, by = "V1", sort = F, all.x = T)  # note order will change
  } else {
    expr_allGenes <- data.frame(geneID, V2 = 0, stringsAsFactors = F)
  }
  expr_value <- expr_allGenes[, "V2"]
  names(expr_value) <- expr_allGenes[, "V1"]
  expr_value <- expr_value[geneID$V1]
  return(expr_value)
})
stopCluster(cl); rm(cl)
expr_rmdup_ensmebl[is.na(expr_rmdup_ensmebl)] <- 0
write.table(x = expr_rmdup_ensmebl, file = paste0(OUT, "/", batch_id, "/exprRmdup/UMIcount_ensemblGene.txt"), 
            row.names = T, col.names = NA, quote = F, sep = "\t")
# ensembl to symbol
gID <- read.table("../../Genomes/human/gene_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)
eID <- read.table("../../Genomes/human/ERCC_ID2Name.txt", header = F, sep = "\t", stringsAsFactors = F)
xID <- rbind(gID, eID)
colnames(xID) <- c("ensembl", "symbol")
expr_rmdup_symbol <- merge(xID, expr_rmdup_ensmebl, by.x = "ensembl", by.y = 0, sort = F)
# remove dup symbol
dupName <- read.table("../../Genomes/human/gene_dupName.txt", header = F , sep = "\t", stringsAsFactors = F)[, 1]
expr_rmdup_symbol <- subset(expr_rmdup_symbol, ! (symbol %in% dupName))
expr_rmdup_symbol <- data.frame(expr_rmdup_symbol[, - c(1, 2)], row.names = expr_rmdup_symbol$symbol, check.names = F)
expr_rmdup_symbol_gene <- expr_rmdup_symbol[rownames(expr_rmdup_symbol) %in% gID$V2, ]
expr_rmdup_symbol_ercc <- expr_rmdup_symbol[rownames(expr_rmdup_symbol) %in% eID$V2, ]
write.table(x = expr_rmdup_symbol_gene, file = paste0(OUT, "/", batch_id, "/exprRmdup/UMIcount_allGenes.txt"), 
            row.names = T, col.names = NA, quote = F, sep = "\t")
write.table(x = expr_rmdup_symbol_ercc, file = paste0(OUT, "/", batch_id, "/exprRmdup/UMIcount_ERCC.txt"), 
            row.names = T, col.names = NA, quote = F, sep = "\t")

### 5. conta_res to stat ----
deConta_stat_LS <- lapply(deConta_res, "[[", 2)
deConta_stat_LS <- lapply(seq_along(deConta_stat_LS), function(i) {
  y <- deConta_stat_LS[[i]]
  y$gene <- names(deConta_stat_LS)[i]
  return(y)
})
deConta_stat_DF <- do.call("rbind", deConta_stat_LS)
deConta_stat_DF$cell <- factor(deConta_stat_DF$cell, levels = tab_DF_sub$cell)
rownames(deConta_stat_DF) <- NULL
deConta_stat_agg <- aggregate(x = deConta_stat_DF[, 2:11], by = list(cell = deConta_stat_DF$cell), FUN = sum)
# molecule level
deConta_stat_agg$ml_ratioRow <- deConta_stat_agg$ml_contaRow / deConta_stat_agg$ml_total
deConta_stat_agg$ml_ratioCol <- deConta_stat_agg$ml_contaCol / deConta_stat_agg$ml_total
deConta_stat_agg$ml_ratioOth <- deConta_stat_agg$ml_contaOth / deConta_stat_agg$ml_total
# type level
deConta_stat_agg$tp_ratioRow <- deConta_stat_agg$tp_contaRow / deConta_stat_agg$tp_total
deConta_stat_agg$tp_ratioCol <- deConta_stat_agg$tp_contaCol / deConta_stat_agg$tp_total
deConta_stat_agg$tp_ratioOth <- deConta_stat_agg$tp_contaOth / deConta_stat_agg$tp_total

write.table(x = deConta_stat_agg, file = paste0(OUT, "/", batch_id, "/", "interCellConta_stat.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")
