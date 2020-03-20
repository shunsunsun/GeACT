# pool datasets from different seqID
setwd("~/lustre/06-Human_cell_atlas/pooled_data/13_heart/")

library("parallel")

# options
dts <- data.frame(seqid = c("21_200110_50","21_200110_50","21_200110_50","22_200120_46","22_200120_46","22_200120_46"), 
                  posid = c("BM_D","FDM_D","ZDM_D","BM_D","FDM_D","ZDM_D"), stringsAsFactors = F)
#out <- "X"

# 1. UMI all genes
cl <- makeCluster(min(20, nrow(dts)))
dts_LS <- parLapply(cl, split(dts, 1:nrow(dts)), function(x) { 
  print(x)
  y0 <- read.table(paste0("../../datasets/", x[1], "/03-expression/merged/UMIcount_allGenes.txt"), header = T, sep = "\t", 
                  stringsAsFactors = F, row.names = 1, check.names = F, comment.char = "")
  y1 <- y0[, grep(paste0("^", x[2], "-"), colnames(y0))]
  return(y1)
})
stopCluster(cl); rm(cl)
names(dts_LS) <- NULL
dts_DF <- do.call("cbind", dts_LS)
dir.create(path = paste0("03-expression/merged"), showWarnings = F, recursive = T)
data.table::fwrite(x = dts_DF, file = paste0("03-expression/merged/UMIcount_allGenes.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)

# 2. cleanFq stat
mt1_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y0 <- read.table(paste0("../../datasets/", x[1], "/01-cleandata/merged/cleanFqStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y1 <- y0[grep(paste0("^", x[2], "-"), y0[, 5]), ]
  y1$seqid <- as.character(x[1])
  return(y1)
})
names(mt1_LS) <- NULL
mt1_DF <- do.call("rbind", mt1_LS)
dir.create(path = paste0("01-cleandata/merged"), showWarnings = F, recursive = T)
write.table(x = mt1_DF, file = paste0("01-cleandata/merged/cleanFqStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 3. map stat
mt2_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y0 <- read.table(paste0("../../datasets/", x[1], "/02-alignment/merged/mapStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y1 <- y0[grep(paste0("^", x[2], "-"), y0[, 2]), ]
  return(y1)
})
names(mt2_LS) <- NULL
mt2_DF <- do.call("rbind", mt2_LS)
dir.create(path = paste0("02-alignment/merged"), showWarnings = F, recursive = T)
write.table(x = mt2_DF, file = paste0("02-alignment/merged/mapStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

mt3_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y0 <- read.table(paste0("../../datasets/", x[1], "/02-alignment/merged/readsDistri.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y1 <- y0[grep(paste0("^", x[2], "-"), y0[, 2]), ]
  return(y1)
})
names(mt3_LS) <- NULL
mt3_DF <- do.call("rbind", mt3_LS)
dir.create(path = paste0("02-alignment/merged"), showWarnings = F, recursive = T)
write.table(x = mt3_DF, file = paste0("02-alignment/merged/readsDistri.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

# 4. expr stat
mt4_LS <- lapply(split(dts, 1:nrow(dts)), function(x) {
  #print(x)
  y0 <- read.table(paste0("../../datasets/", x[1], "/03-expression/merged/exprStat.txt"), header = F, sep = "\t", stringsAsFactors = F)
  y1 <- y0[grep(paste0("^", x[2], "-"), y0[, 2]), ]
  return(y1)
})
names(mt4_LS) <- NULL
mt4_DF <- do.call("rbind", mt4_LS)
dir.create(path = paste0("03-expression/merged"), showWarnings = F, recursive = T)
write.table(x = mt4_DF, file = paste0("03-expression/merged/exprStat.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
