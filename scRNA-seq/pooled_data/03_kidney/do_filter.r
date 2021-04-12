# filter cells and genes
setwd("~/lustre/06-Human_cell_atlas/pooled_data/03_kidney/")

library("reshape2")
library("ggplot2")
library("cowplot")
source("../../scripts/filter_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/filtering/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/filtering.RData"))

# load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Genomes/human/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(geneID)
colnames(geneID) <- c("ensembl", "symbol")
geneID$ensembl_alt <- gsub("\\.[0-9]+", "", geneID$ensembl)

# load expression data
exprMatr <- read.table("03-expression/merged/UMIcount_allGenes.txt",header = T,row.names = 1,sep = "\t", stringsAsFactors = F, check.names = F, comment.char = "")
dim(exprMatr)

# extract data
#exprMatr <- exprMatr[, grep(paste0("^", samplingPos, "-"), colnames(exprMatr))]

# remove control/abnormal sample
black_list <- "nothing"
exprMatr <- exprMatr[, grep(black_list, colnames(exprMatr), invert = T)]

# remove 4th inner
#inner_id <- as.numeric(gsub(".*_", "", colnames(exprMatr))) %% 8; inner_id[inner_id==0] <- 8; inner_id <- factor(inner_id)
#exprMatr <- exprMatr[, inner_id!=4]

print("Before filter:")
dim(exprMatr)

## 1. cell level ----
# 1 clean reads and A/B ratio
cleanStat_per_cell <- read.table("01-cleandata/merged/cleanFqStat.txt", header = F, sep = "\t", row.names = 5, stringsAsFactors = F, check.names = F)
dim(cleanStat_per_cell)
cleanStat_per_cell <- merge(x = colnames(exprMatr), y = cleanStat_per_cell, by.x = 1, by.y = 0, sort = F)
rownames(cleanStat_per_cell) <- cleanStat_per_cell$x
cleanStat_per_cell <- cleanStat_per_cell[, -1]
cleanStat_per_cell$ABratio <- cleanStat_per_cell$V8/(cleanStat_per_cell$V8 + cleanStat_per_cell$V9)
hist(cleanStat_per_cell$V12/1e6, breaks = 40, xlab = "Clean reads (M)", main = NA)
hist(cleanStat_per_cell$ABratio, breaks = 40, xlab = "A/B ratio", main = NA)
cellStat <- cleanStat_per_cell[, c("ABratio", "V12")]
colnames(cellStat)[2] <- "cleanReads"
# 2 clean reads and mapping ratio
mapStat_per_cell <- read.table("02-alignment/merged/mapStat.txt",header = F, sep = "\t", row.names = 2, stringsAsFactors = F, check.names = F)
dim(mapStat_per_cell)
mapStat_per_cell <- merge(x = colnames(exprMatr), y = mapStat_per_cell, by.x = 1, by.y = 0, sort = F)
rownames(mapStat_per_cell) <- mapStat_per_cell$x
mapStat_per_cell <- mapStat_per_cell[, -1]
hist(mapStat_per_cell$V7,breaks = 40, xlab = "Mapping ratio", main = NA, xlim = c(0, 1))
cellStat$mpRatio <- mapStat_per_cell$V7
# 3 detected genes
expressed_genes_per_cell <- apply(exprMatr>0, 2, sum)
hist(expressed_genes_per_cell,breaks = 40, xlab = "Detected genes", main = NA)
cellStat$nGene <- expressed_genes_per_cell
# 4 UMI
nUMI_per_cell <- colSums(exprMatr)
hist(nUMI_per_cell/1e3,breaks = 40, xlab = "Detected transcripts (k)", main = NA)
cellStat$nUMI <- nUMI_per_cell
# 5 mito ratio
#mito.genes <- grep(pattern = "^MT-", x = rownames(exprMatr), value = T)
percent.mito <- colSums(exprMatr[grep(pattern = "^MT-", x = rownames(exprMatr)),])/colSums(exprMatr)
hist(percent.mito, breaks = 40, xlab = "Mito. ratio", main = NA)
cellStat$mitoRatio <- percent.mito
#save.image(file="filtering.RData")
#load("filtering.RData")
# 6 ERCC ratio
exprStat_per_cell <- read.table("03-expression/merged/exprStat.txt", header = F, sep = "\t", row.names = 2, stringsAsFactors = F)
dim(exprStat_per_cell)
exprStat_per_cell <- merge(x = colnames(exprMatr), y = exprStat_per_cell, by.x = 1, by.y = 0, sort = F)
hist(exprStat_per_cell$V6 / exprStat_per_cell$V3, breaks = 40, xlab = "ERCC ratio", main = NA)
cellStat$ERCCratio <- exprStat_per_cell$V6 / exprStat_per_cell$V3
# 7 correlation
# cor_cells <- cor(exprMatr[rowSums(exprMatr)>0, ], method = "spearman")
# diag(cor_cells) <- NA
# avg_cor <- rowMeans(cor_cells, na.rm = T)
# hist(avg_cor, breaks = 40, xlab = "Pairwise correlation", main = NA)
cellStat$avgCor <- 1
# 8 doublets
library("mgcv")
plot(cellStat$nUMI, cellStat$nGene, xlab = "UMI number", ylab = "Gene number")
gam.reg <- gam(nUMI ~ s(nGene, bs = "cs"), data = cellStat)
gam.pre <- predict(gam.reg, newdata = list(nGene=cellStat$nGene))
gam_DF <- data.frame(cellStat, pred_value = gam.pre, obs_fc = cellStat$nUMI / gam.pre)
cellStat$doublet <- gam_DF$nUMI>quantile(gam_DF$nUMI, 0.99) & gam_DF$obs_fc>2
table(cellStat$doublet)

# merged
do_plotIndex()

plot(cellStat$nUMI, cellStat$nGene)
plot(cellStat$cleanReads/1e6, cellStat$ABratio, xlab = "Clean reads (M)", ylab = "A/B ratio"); abline(v = 0.5, col = "blue", lty = 2)
plot(cellStat$cleanReads/1e6, cellStat$mitoRatio, xlab = "Clean reads (M)", ylab = "Mito. ratio"); abline(v = 0.5, col = "blue", lty = 2)
plot(cellStat$cleanReads/1e6, cellStat$ERCCratio, xlab = "Clean reads (M)", ylab = "ERCC ratio"); abline(v = 0.5, col = "blue", lty = 2)
plot(cellStat$cleanReads/1e6, cellStat$avgCor, xlab = "Clean reads (M)", ylab = "Pairwise correlation"); abline(v = 0.5, col = "blue", lty = 2)

### setting cutoff for filtering
cutoff_DF <- data.frame(ABratio = 0.9, cleanReads = 0.4*1e6, mpRatio = 0.6, 
                        nGene_l = 1000, nGene_h = 10000, nUMI_l = 3*1e3, nUMI_h = 124*1e3,
                        mitoRatio = 0.2, ERCCratio = 0.25, avgCor = 0.15)
###

filtering_out <- do_cellFiltering()
# more filtering
# filtering_out$res[cellStat$cleanReads > 20e6] <- F

cellStat$filter <- filtering_out$res
# report
do_show_ftd_stat(cellStat)

sink(paste0(OUT, "/filtering_stat.txt"))
print(cutoff_DF)
cat("\n")
do_show_ftd_stat(cellStat)
sink()

# plot
pdf(paste0(OUT, "/filtering_cells.pdf"), width = 4, height = 4, useDingbats = F)
library("ggplot2")

p <- ggplot(filtering_out$tb, aes(x = cell, y = index, fill = as.factor(value))) + geom_tile(show.legend = F) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(size = 1, color = "black")) + 
  scale_y_discrete(limits = rev(levels(filtering_out$tb$index)), position = "right", expand = c(0, 0))
print(p)

# stat for each bigBatch
cellStat_DF <- cellStat
cellStat_DF$bigBatch <- gsub("-.*", "", rownames(cellStat_DF))
cellStat_DF_melted <- melt(table(cellStat_DF[, c("bigBatch", "filter")]))
ggplot(cellStat_DF_melted, aes(x = bigBatch, y = value, fill = bigBatch, alpha = filter)) + 
  geom_bar(stat = "identity", show.legend = F) + 
  scale_alpha_discrete(range = c(0.4,1)) + 
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(0, max(aggregate(cellStat_DF_melted$value, by = list(cellStat_DF_melted$bigBatch), sum)[,2])*1.05), 
                     expand = c(0, 0)) + 
  ylab("Cell number") + 
  geom_text(aes(label = value, y = value), size = 3, position = position_stack(vjust = 0.5), show.legend = F)

# stat for each feature
ggplot(cellStat, aes(ABratio, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("A/B ratio") + ylab("Cell number") + 
  geom_vline(xintercept = cutoff_DF$ABratio, linetype = "dashed")

ggplot(cellStat, aes(cleanReads/1e6, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("Clean reads number / (Million)") + ylab("Cell number") + 
  geom_vline(xintercept = cutoff_DF$cleanReads/1e6, linetype = "dashed")

ggplot(cellStat, aes(mpRatio, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("Mapping ratio") + ylab("Cell number") + 
  geom_vline(xintercept = cutoff_DF$mpRatio, linetype = "dashed")

ggplot(cellStat, aes(nGene, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("Dected gene number") + ylab("Cell number") + 
  geom_vline(xintercept = c(cutoff_DF$nGene_l, cutoff_DF$nGene_h), linetype = "dashed")

tmp <- data.frame(cellStat, bigBatch = gsub("-.*", "", rownames(cellStat)), stringsAsFactors = F)
if(length(unique(tmp$bigBatch))>1) {
ggplot(tmp, aes(nGene, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black", show.legend = F) + 
  xlab("Dected gene number") + ylab("Cell number") + facet_grid(bigBatch ~ .) + 
  geom_vline(xintercept = c(cutoff_DF$nGene_l, cutoff_DF$nGene_h), linetype = "dashed")
}

ggplot(cellStat, aes(nUMI, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("UMI count") + ylab("Cell number") + 
  geom_vline(xintercept = c(cutoff_DF$nUMI_l, cutoff_DF$nUMI_h), linetype = "dashed")

ggplot(cellStat, aes(log10(nUMI), fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab(expression(paste(log[10], " (UMI count)"))) + ylab("Cell number") + 
  geom_vline(xintercept = log10(c(cutoff_DF$nUMI_l, cutoff_DF$nUMI_h)), linetype = "dashed")

ggplot(cellStat, aes(mitoRatio, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("Expression raio of mitochondrial genes") + ylab("Cell number") + 
  geom_vline(xintercept = cutoff_DF$mitoRatio, linetype = "dashed")

ggplot(cellStat, aes(avgCor, fill=filter)) + geom_histogram(bins = 40, alpha=.5, position="identity", color="black") + 
  theme_grey() + xlab("Average of pairwise Spearman correlation") + ylab("Cell number") + 
  geom_vline(xintercept = cutoff_DF$avgCor, linetype = "dashed")

ggplot(cellStat, aes(x = cleanReads/1e6, y = ABratio, color=filter, size=nUMI)) + geom_point(alpha=0.6) + 
  theme_grey() + xlab("Reads number / (Million)") + ylab("A/B ratio") + 
  geom_vline(xintercept = cutoff_DF$cleanReads/1e6, linetype = "dashed") + 
  geom_hline(yintercept = cutoff_DF$ABratio, linetype = "dashed")

ggplot(cellStat, aes(x = cleanReads/1e6, y = ERCCratio, color=filter, size=nUMI)) + geom_point(alpha=0.6) + 
  theme_grey() + xlab("Reads number / (Million)") + ylab("ERCC ratio") + 
  geom_vline(xintercept = cutoff_DF$cleanReads/1e6, linetype = "dashed") + 
  geom_hline(yintercept = cutoff_DF$ERCCratio, linetype = "dashed")

ggplot(cellStat, aes(x = cleanReads/1e6, y = avgCor, color=filter, size=nUMI)) + geom_point(alpha=0.6) + 
  theme_grey() + xlab("Reads number / (Million)") + ylab("Average Spearman correlation") + 
  geom_vline(xintercept = cutoff_DF$cleanReads/1e6, linetype = "dashed") + 
  geom_hline(yintercept = cutoff_DF$avgCor, linetype = "dashed")

ggplot(cellStat, aes(x = nUMI, y = nGene, color=filter, size=cleanReads, shape=doublet==1)) + geom_point(alpha=0.6) + 
  theme_grey() + xlab("UMI number") + ylab("Detected gene number") + 
  geom_vline(xintercept = c(cutoff_DF$nUMI_l, cutoff_DF$nUMI_h), linetype = "dashed") + 
  geom_hline(yintercept = c(cutoff_DF$nGene_l, cutoff_DF$nGene_h), linetype = "dashed") + 
  scale_shape_discrete(name="doublet") + 
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"), color = "1", show.legend = F, alpha = 0.6)

dev.off()

## 2. gene level ----
# 1
expressed_cells_per_gene <- apply(exprMatr[, cellStat$filter]>0, 1, sum)
hist(expressed_cells_per_gene/ncol(exprMatr[, cellStat$filter]), breaks = 40, xlab = "Ratio of cells expressed", main = NA)
abline(v=0.05,lty=2,col="blue")
# 2
nCountsPerGene <- apply(exprMatr, 1, sum)
hist(nCountsPerGene, breaks = 1000000, xlim = c(0,1000))
# merge
geneStat <- data.frame(nCell=expressed_cells_per_gene, nUMI=nCountsPerGene)
geneStat$validCell <- sum(cellStat$filter)
geneStat$cellRatio <- geneStat$nCell/geneStat$validCell
geneStat$filter <- geneStat$nCell>=10
table(geneStat$filter)

## 3. filtering both of them
exprMatr_filtered <- exprMatr[geneStat$filter, cellStat$filter]
print("After filter:")
dim(exprMatr_filtered)

# filtering only cells
exprMatr_cellFiltered <- exprMatr[, cellStat$filter]
dim(exprMatr_cellFiltered)

exprMatr_cellFiltered_CPM <- sweep(exprMatr_cellFiltered, MARGIN = 2, STATS = colSums(exprMatr_cellFiltered), FUN = "/") * 1e6
exprMatr_cellFiltered_CPM_ensembl <- exprMatr_cellFiltered_CPM
rownames(exprMatr_cellFiltered_CPM_ensembl) <- geneID$ensembl_alt[match(rownames(exprMatr_cellFiltered_CPM_ensembl), geneID$symbol)]

# output
data.table::fwrite(x = exprMatr, file = paste0(OUT, "/UMIcount_unfiltered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", paste0(OUT, "/UMIcount_unfiltered.txt"), " > ", OUT, "/UMIcount_unfiltered.txt.gz"))

data.table::fwrite(x = exprMatr_filtered, file = paste0(OUT, "/UMIcount_filtered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", paste0(OUT, "/UMIcount_filtered.txt"), " > ", OUT, "/UMIcount_filtered.txt.gz"))

data.table::fwrite(x = exprMatr_cellFiltered, file = paste0(OUT, "/UMIcount_cellFiltered.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
system(paste0("gzip -c ", paste0(OUT, "/UMIcount_cellFiltered.txt"), " > ", OUT, "/UMIcount_cellFiltered.txt.gz"))

data.table::fwrite(x = exprMatr_cellFiltered_CPM, file = paste0(OUT, "/UMIcount_cellFiltered_CPM.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)
data.table::fwrite(x = exprMatr_cellFiltered_CPM_ensembl, file = paste0(OUT, "/UMIcount_cellFiltered_CPM_ensembl.txt"), row.names = T, col.names = T, quote = F, sep = "\t", nThread = 10)

write.table(x = cellStat,file = paste0(OUT, "/filtering_cells.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")
write.table(x = geneStat,file = paste0(OUT, "/filtering_genes.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")

# extract normalized data from filtered data
# exprMatr_CPM <- sweep(x = exprMatr, MARGIN = 2, STATS = colSums(exprMatr), FUN = "/") * 0.1 * 1e6 # scale by 0.1M
# exprMatr_normed <- exprMatr_CPM[rownames(exprMatr_filtered), cellStat$filter]
# write.table(x = exprMatr_normed, file = paste0(OUT, "/UMIcount_normed.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")

# log2 transformation
# exprMatr_normed_log2 <- log2(exprMatr_normed + 1)
# write.table(x = exprMatr_normed_log2, file = paste0(OUT, "/UMIcount_normed_log2.txt"), row.names = T, col.names = NA, quote = F, sep = "\t")

save.image(file = paste0(OUT, "/filtering.RData"))
