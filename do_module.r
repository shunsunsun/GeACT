#!/usr/bin/Rscript

library("GENIE3")

exprMatr <- read.table(file = "04-expression/merged/UMIcount_filtered.txt",header = T,row.names = 1,sep = "\t")
rownames(exprMatr) <- sapply(strsplit(x = rownames(exprMatr), split = ".", fixed = T),"[",1)
exprMatr <- as.matrix(exprMatr)
dim(exprMatr)

TF_list <- read.table(file = "Data/human/TF_list.txt", header = F,stringsAsFactors = F)[,1]
TF_in_table <- TF_list[TF_list%in%rownames(exprMatr)]

set.seed(123)	# For reproducibility of results
weightMat <- GENIE3(exprMatrix = exprMatr, regulators = TF_in_table, nCores = 20, verbose = TRUE)	# slow
dim(weightMat)

linkList <- getLinkList(weightMat, threshold = 0.001)
dim(linkList)
head(linkList)

