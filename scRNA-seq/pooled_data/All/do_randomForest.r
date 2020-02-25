# random forest analysis
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("ggplot2")
library("cowplot")
library("randomForest")
#source("../../scripts/cluster_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/randomForest/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/randomForest.RData"))

TF <- read.table("../../Data/human_TF.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(TF)
lncRNA <- read.table("../../Data/human_lncRNA.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(lncRNA)
chromtainRM <- read.table("../../Data/human_cr.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(chromtainRM)

# 1. pre-process ----
# Load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Genomes/human/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(geneID)
colnames(geneID) <- c("ensembl", "symbol")

# Load the expr matrix
expr_data <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F, comment.char = "")
dim(expr_data)
# norm
expr_data_normed <- sweep(expr_data, 2, colSums(expr_data), "/")

# Load the cell metatable
cellMetaData <- read.table(file = "cell_metatable_filtered.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
dim(cellMetaData)
all(colnames(expr_data) == rownames(cellMetaData))
cellMetaData$ts_ident <- Hmisc::capitalize(paste(cellMetaData$tissue, cellMetaData$ident, sep = "."))
#cellMetaData$ts_clgrp <- Hmisc::capitalize(paste(cellMetaData$tissue, ident2clgrp(cellMetaData$ident), sep = "."))
#length(unique(cellMetaData$ts_clgrp))

# avg
#expr_data_avg <- sapply(split(rownames(cellMetaData), cellMetaData$ts_ident), function(x) { y <- rowMeans(expr_data_normed[, x, drop = F]) })

# extract TF
expr_data_normed_onlyTF <- as.data.frame(t(expr_data_normed[rownames(expr_data_normed) %in% TF$V1, ]))
# modify special char
colnames(expr_data_normed_onlyTF) <- gsub("-", "___", colnames(expr_data_normed_onlyTF))
# add cell type
expr_data_normed_onlyTF$cellType <- cellMetaData$ts_ident

## all cell types ----
# split data
set.seed(1234)
split_idx <- sample(nrow(expr_data_normed_onlyTF), nrow(expr_data_normed_onlyTF) * 0.7)
expr_data_train <- expr_data_normed_onlyTF[split_idx, ]
expr_data_train$cellType <- factor(expr_data_train$cellType)
expr_data_test <- expr_data_normed_onlyTF[- split_idx, ]

set.seed(1234)
#rf_fit <- randomForest(cellType ~ ., data = expr_data_train, importance = T)
varImpPlot(rf_fit)

rf_importance <- importance(rf_fit, type = 2)
rf_importance <- as.data.frame(rf_importance)
rownames(rf_importance) <- gsub("___", "-", rownames(rf_importance))
rf_importance$gene <- rownames(rf_importance)
rf_importance <- rf_importance[order(rf_importance$MeanDecreaseGini, decreasing = T), ]

rf_importance_sub <- head(rf_importance, 10)
rf_importance_sub$gene <- factor(rf_importance_sub$gene, levels = rev(unique(rf_importance_sub$gene)))

pdf(paste0(OUT, "/rf_TF_allCellType.pdf"), width = 5, height = 3.5)

ggplot(rf_importance_sub, aes(x = gene, y = MeanDecreaseGini)) + 
  geom_point(color = NA) + 
  geom_vline(xintercept = 1:nrow(rf_importance_sub), linetype = "dashed", color = "grey") + 
  geom_point(color = "dodgerblue", size = 3) + 
  xlab("TF") + ylab("Mean decrease Gini") + 
  #ggtitle("All cell types") + 
  coord_flip()

dev.off()

## Epi ----
cells_sub <- rownames(cellMetaData)[grepl("^Epi", cellMetaData$ident) | (cellMetaData$ident %in% c("PT", "LoH", "LoH-Prog", "DT", "PC-CLU", "PC-BCAT1", "Podocyte-GPC3", "Podocyte-PLA2R1"))]
expr_data_normed_onlyTF <- expr_data_normed_onlyTF[cells_sub, ]
table(expr_data_normed_onlyTF$cellType)

# split data
set.seed(1234)
split_idx <- sample(nrow(expr_data_normed_onlyTF), nrow(expr_data_normed_onlyTF) * 0.7)
expr_data_train <- expr_data_normed_onlyTF[split_idx, ]
expr_data_train$cellType <- factor(expr_data_train$cellType)
expr_data_test <- expr_data_normed_onlyTF[- split_idx, ]

set.seed(1234)
#rf_fit2 <- randomForest(cellType ~ ., data = expr_data_train, importance = T)
varImpPlot(rf_fit2)

rf_importance <- importance(rf_fit2, type = 2)
rf_importance <- as.data.frame(rf_importance)
rownames(rf_importance) <- gsub("___", "-", rownames(rf_importance))
rf_importance$gene <- rownames(rf_importance)
rf_importance <- rf_importance[order(rf_importance$MeanDecreaseGini, decreasing = T), ]

rf_importance_sub <- head(rf_importance, 10)
rf_importance_sub$gene <- factor(rf_importance_sub$gene, levels = rev(unique(rf_importance_sub$gene)))

pdf(paste0(OUT, "/rf_TF_Epi.pdf"), width = 5, height = 3)

ggplot(rf_importance_sub, aes(x = gene, y = MeanDecreaseGini)) + 
  geom_point(color = NA) + 
  geom_vline(xintercept = 1:nrow(rf_importance_sub), linetype = "dashed", color = "grey") + 
  geom_point(color = "dodgerblue", size = 3) + 
  xlab("TF") + ylab("Mean decrease Gini") + 
  ggtitle("Epithelium") + 
  coord_flip()

dev.off()

save.image(file = paste0(OUT, "/randomForest.RData"))
