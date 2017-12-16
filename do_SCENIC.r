# GENIE3 ----
library("GENIE3")
library("reshape2")

setwd("/home/tianf/lustre/06-Human_cell_atlas/test_samples")

# input ----
exprMatr <- read.table(file = "03-expression/merged/UMIcount_filtered.txt",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
exprMatr <- as.matrix(exprMatr)
dim(exprMatr)
# within sample normalization
exprMatr_CPM <- sweep(x = exprMatr, MARGIN = 2, STATS = colSums(exprMatr), FUN = "/") * 1e6
exprMatr_logCPM <- log2(exprMatr_CPM + 1)

# TF from Animal TFBS
gene_ID2name <- read.table(file = "Genomes/human/gene_ID2Name.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors = F)
rownames(gene_ID2name) <- gsub("\\.[0-9]*","",rownames(gene_ID2name))
TF_list_ensmebl <- read.table(file = "Data/human/TF_list.txt", header = F,stringsAsFactors = F)[,1]
length(TF_list_ensmebl)
TF_list_ensmebl_in_table <- TF_list_ensmebl[TF_list_ensmebl%in%rownames(gene_ID2name)]  # TF in gencode annotation
length(TF_list_ensmebl_in_table)

TF_list <- gene_ID2name[TF_list_ensmebl_in_table,]
#write.table(file = "03-expression/merged/TF_list_id.txt", x = TF_list, row.names = F, col.names = F, quote = F, sep = "\t")
TF_list_in_table <- TF_list[TF_list%in%rownames(exprMatr)]  # TF in expression table
length(TF_list_in_table)

# TF from Rcistarget
library(RcisTarget.hg19.motifDatabases.20k)
data(hg19_direct_motifAnnotation)
allTFs <- hg19_direct_motifAnnotation$allTFs
length(allTFs)
allTFs <- allTFs[allTFs%in%rownames(exprMatr)]
length(allTFs)

# GENIE3 ----
set.seed(123) # For reproducibility of results
library("foreach")
#weightMat <- GENIE3(exprMatrix = exprMatr, regulators = allTFs, nCores = 50, verbose = TRUE) # slow
save(weightMat, file="03-expression/merged/GENIE3_weightMatrix.RData")
dim(weightMat)

linkList <- getLinkList(weightMat, threshold = 0.001)
dim(linkList)
colnames(linkList) <- c("TF", "Target", "weight")
head(linkList)
save(linkList, file="03-expression/merged/GENIE3_linkList.RData")
#linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),] # order by weight

# Creating TF modules (potential TF-targets)
quantile(linkList$weight, probs=c(0.75, 0.90))

plot(linkList$weight[1:2000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
     ylab="Weight", xlab="Links sorted decreasingly")
#abline(h=0.001, col="blue") # Threshold
abline(h=0.0015, col="blue")

sum(linkList$weight>0.0015)/nrow(linkList)
linkList_cut <- linkList[which(linkList[,"weight"]>0.0015),]
# Number of links over the threshold: 
nrow(linkList_cut) 

linkList_cut$TF <- as.character(linkList_cut$TF)
linkList_cut$Target <- as.character(linkList_cut$Target)

# Create the gene-sets / TF-modules:
tfModules <- linkList_cut
rbind(nGeneSets=nrow(tfModules), nTFs=length(unique(tfModules$TF)), nTargets=length(unique(tfModules$Target)))
save(tfModules, file="03-expression/merged/GENIE3_tfModules.RData")

# correlation
source("bigcor.r")
corrMat_tmp <- bigcor(x = t(exprMatr), nblocks = 5, mt = "spearman", cpu.num = 20, do.melt = T, do.merge = T)
corrMat <- dcast(data = corrMat_tmp, formula = Var1 ~ Var2)
rownames(corrMat) <- corrMat$Var1
corrMat <- corrMat[, -1]
dim(corrMat)
save(corrMat, file="03-expression/merged/GENIE3_corrMat.RData")

# correlation after normalization
corrMat_tmp_lc <- bigcor(x = t(exprMatr_logCPM), nblocks = 5, mt = "spearman", cpu.num = 20, do.melt = T, do.merge = T)
corrMat_lc <- dcast(data = corrMat_tmp_lc, formula = Var1 ~ Var2)
rownames(corrMat_lc) <- corrMat_lc$Var1
corrMat_lc <- corrMat_lc[, -1]
dim(corrMat_lc)
save(corrMat_lc, file="03-expression/merged/GENIE3_corrMat_lc.RData")

hist(as.matrix(corrMat)[,1:100], breaks = 100)
hist(as.matrix(corrMat_lc)[,1:100], breaks = 100)

# partial correlation
#library("ppcor")
#expr_pcor_LS <- pcor(x = t(exprMatr), method = "spearman")
library("psych")
partial_cor <- partial.r(corrMat, c("E2F7", "E2F1"), which(! colnames(corrMat)%in%c("E2F7", "E2F1")), method = "spearman")

# Split into positive- and negative-correlated targets
# Keep only correlation between TFs and potential targets
tfs <- unique(tfModules$TF)
corrMat <- corrMat[tfs,]
dim(corrMat)

# Split TF modules according to correlation
tfModules_byTF <- split(tfModules, factor(tfModules$TF))
library("data.table")
tfModules_withCorr_byTF <- lapply(tfModules_byTF, function(tfGeneSets)
{
  tf <- unique(tfGeneSets$TF)
  targets <- tfGeneSets$Target
  cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
})
tfModules_withCorr <- data.frame(rbindlist(tfModules_withCorr_byTF))
save(tfModules_withCorr, file="03-expression/merged/GENIE3_tfModules_withCorr.RData")


# Rcistarget ----
org <- "human"
if(org=="human")
{
  library(RcisTarget.hg19.motifDatabases.20k)
  
  # Motif rankings (genes x motifs)
  data(hg19_500bpUpstream_motifRanking)
  data(hg19_10kbpAroundTss_motifRanking)
  motifRankings <- list()
  motifRankings[["500bp"]] <- hg19_500bpUpstream_motifRanking
  #motifRankings[["10kbp"]] <- hg19_10kbpAroundTss_motifRanking
  
  # Motif annotation (TFs)
  data(hg19_direct_motifAnnotation)
  direct_motifAnnotation <- hg19_direct_motifAnnotation
  data(hg19_inferred_motifAnnotation) # optional
  inferred_motifAnnotation <- hg19_inferred_motifAnnotation
}

# Remove genes missing from RcisTarget databases
#  (In case the input matrix wasn't already filtered)
tfModules_withCorr <- tfModules_withCorr[which(as.character(tfModules_withCorr$TF) %in% allTFs),]
geneInDb <- tfModules_withCorr$Target %in% motifRankings[["500bp"]]@rankings$rn
# Genes in co-expression modules not available in RcisTargetDatabases:
missingGenes <- sort(unique(tfModules_withCorr[which(!geneInDb),"Target"]))
length(missingGenes)
head(missingGenes)
tfModules_withCorr <- tfModules_withCorr[which(geneInDb),]

# Targets with positive correlation
tfModules_Selected <- tfModules_withCorr[which(tfModules_withCorr$corr==1),]
dim(tfModules_Selected)

# Add a column with the geneSet name (TF_method)
tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, "positive", sep="_"))
head(tfModules_Selected)

# Split into tfModules (TF-modules, with several methods)
tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)

# Keep gene sets with at least 20 genes
length(tfModules)
tfModules <- tfModules[which(lengths(tfModules)>=20)]
length(tfModules)

# Add TF to the gene set (used in the following steps, careful if editing)
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
  tf <- strsplit(gsn, "_")[[1]][1]
  unique(c(tf, tfModules[[gsn]]))
}), names(tfModules))
save(tfModules, file="03-expression/merged/SCENIC_tfModules_forMotifEnrichmet.RData")

tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
sort(table(tfModulesSummary[,2]))

library(RcisTarget)

################################################################
# 1. Calculate motif enrichment for each TF-module

### 1.1 Calculate enrichment
motifs_AUC <- lapply(motifRankings, function(ranking) calcAUC(tfModules, ranking, aucMaxRank=0.01*nrow(ranking@rankings), nCores=30, verbose=FALSE))
save(motifs_AUC, file="03-expression/merged/SCENIC_motifs_AUC_500bp_10kbp.RData")

### 1.2 Conver to table, filter by NES & add the TFs to which the motif is annotated
# (For each database...)
motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
{
  # Extract the TF of the gene-set name (i.e. MITF_w001):
  tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
  
  # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
  addMotifAnnotation(aucOutput, highlightTFs=tf, nesThreshold=3.0, digits=3,
                     motifAnnot_direct=direct_motifAnnotation,
                     motifAnnot_inferred=inferred_motifAnnotation)
})

# Merge both tables, adding a column that contains the 'motifDb' 
motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
  cbind(motifDb=dbName, motifEnrichment[[dbName]])
}))
save(motifEnrichment, file="03-expression/merged/SCENIC_motifEnrichment.RData")
cat("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))

### 1.3 Keep only the motifs annotated to the initial TF
motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
save(motifEnrichment_selfMotifs, file="03-expression/merged/SCENIC_motifEnrichment_selfMotifs.RData")
cat("Number of motifs annotated to the initial TF: ", nrow(motifEnrichment_selfMotifs))
#rm(motifEnrichment)

################################################################
# 2. Prune targets

motifEnrichment_selfMotifs_wGenes <- lapply(names(motifRankings), function(motifDbName){
  addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifDb==motifDbName],
                      geneSets=tfModules,
                      rankings=motifRankings[[motifDbName]],
                      maxRank=5000, method="aprox", nCores=50)
})

library(data.table)
motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
dim(motifEnrichment_selfMotifs_wGenes)
save(motifEnrichment_selfMotifs_wGenes, file="03-expression/merged/SCENIC_motifEnrichment_selfMotifs_wGenes.RData")

# Save as text:
write.table(motifEnrichment_selfMotifs_wGenes, file="03-expression/merged/SCENIC_MotifEnrichment.tsv", sep="\t", quote=FALSE, row.names=FALSE)

motifEnrichment_selfMotifs_wGenes[order(NES,decreasing=TRUE)][1:5,-"enrichedGenes", with=F]

# Format regulons
library(data.table)
motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
  genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
  oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
  data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
})
motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)

# Get targets for each TF, but keep info about best motif/enrichment 
# (directly annotated motifs are considered better)
regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
  # print(unique(tfTargets$TF))
  tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
    directAnnot <- "**" %in% enrOneGene$annot
    enrOneGeneByAnnot <- enrOneGene
    if(directAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
    bestMotif <- which.max(enrOneGeneByAnnot$NES)
    
    cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene), 
          bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]), 
          directAnnot=directAnnot)
  })), stringsAsFactors=FALSE)
  tfTable[order(tfTable$NES, decreasing = TRUE),]
})
regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "directAnnot")

# Optional: Add Genie3 score
#load("int/1.4_GENIE3_linkList.RData")
#linkList <- linkList[which(linkList$weight>=0.001),]
#rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
#regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])

save(regulonTargetsInfo, file="03-expression/merged/SCENIC_regulonTargetsInfo.RData")
write.table(regulonTargetsInfo, file="03-expression/merged/SCENIC_regulonTargetsInfo.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
# only for the directAnnot
regulonTargetsInfo_directAnnot <- regulonTargetsInfo[regulonTargetsInfo$directAnnot=="TRUE", ]
dim(regulonTargetsInfo_directAnnot)
write.table(regulonTargetsInfo_directAnnot, file="03-expression/merged/SCENIC_regulonTargetsInfo_directAnnot.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Split into regulons (output: list TF –> targets)
regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$directAnnot)
regulons <- sapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
regulons_extended <- sapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(x[,"gene"]))
regulons_extended <- sapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], regulons_extended[[tf]]))))
names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
regulons <- c(regulons, regulons_extended)
save(regulons, file="03-expression/merged/SCENIC_regulons_asGeneSet.RData")

# How many TFs are self-regulating?
table(sapply(names(regulons), function(x) x %in% regulons[[x]]))

# Motifs associated to a TF
#selTF <- "E2F1"
#subsetTable <- motifEnrichment_selfMotifs_wGenes[highlightedTFs %in% selTF][order(NES,decreasing=TRUE)][,-"enrichedGenes", with=F]
#subsetTable <- addLogo(subsetTable)
#library(DT)
#datatable(subsetTable, escape=FALSE, filter="top", options=list(pageLength=5))

# Network activity in each cell
regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
regulons <- regulons[lengths(regulons)>=10]

# Add the TF & rename
regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
save(regulons, file="03-expression/merged/SCENIC_regulons_forAUCell.RData")
length(regulons)
cbind(names(regulons)[1:10])

library(AUCell)
# 1. Create rankings
aucellRankings <- AUCell.buildRankings(exprMatr, nCores=30, plotStats=TRUE)
abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
save(aucellRankings, file="03-expression/merged/SCENIC_aucellRankings.RData")

# 2. Calculate AUC
library("foreach")
regulonAUC <- AUCell.calcAUC(regulons, aucellRankings, aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=30)
save(regulonAUC, file="03-expression/merged/SCENIC_regulonAUC.RData")

# Order the modules by similarity, for easier exploration in the upcoming steps & save
variableRegulons <- names(which(apply(getAuc(regulonAUC), 1, sd) > 0))
reguDist <-as.dist(1-cor(t(getAuc(regulonAUC)[variableRegulons,]), method="spear"))
reguClust <- hclust(reguDist, method="ward.D2")
regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
regulonOrder <- reguClust$labels[reguClust$order]
regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
regulonAUC@matrix <- regulonAUC@matrix[regulonOrder,]
#save(regulonAUC, file="int/3.2_regulonAUC.RData")

# Overview of cell states according to module activity (tSNE on AUC)
# (It is recommended to try different perplexity values)
regulonAUC_subset <- subset(regulonAUC, onlyNonDirectExtended(rownames(regulonAUC)))

### cell type specific active module
cell_type <- read.table("03-expression/merged/seurat_metaData.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
dim(cell_type)
marker_genes <- read.table(file = "03-expression/merged/marker_genes.txt", header = T, sep = "\t", stringsAsFactors = F)

regulonAUC_MT <- regulonAUC_subset@matrix
regulonAUC_MT <- regulonAUC_MT[,rownames(cell_type)]  # sort by cell name in meta table
regulonAUC_cmp <- t(apply(regulonAUC_MT, 1, function(x) {
  y1 <- x[cell_type$ident=="G1/S"]
  y2 <- x[cell_type$ident=="G2/M"]
  pv <- wilcox.test(y1, y2)$p.value
  rt <- mean(y1)/mean(y2)
  z <- c(pv, rt)
}))
regulonAUC_cmp <- as.data.frame(regulonAUC_cmp)
colnames(regulonAUC_cmp) <- c("p-value","ratio")
regulonAUC_cmp <- regulonAUC_cmp[order(regulonAUC_cmp[,"p-value"]), ]
regulonAUC_cmp$`q-value` <- p.adjust(regulonAUC_cmp$`p-value`, method = "BH")
regulonAUC_csp <- regulonAUC_cmp[regulonAUC_cmp$`q-value`<0.05, ]

regulonAUC_sig <- regulonAUC_MT[rownames(regulonAUC_csp)[c(which(regulonAUC_csp$ratio>1), which(regulonAUC_csp$ratio<1))], ]
regulonAUC_sig <- regulonAUC_sig[,c(rownames(cell_type)[cell_type$ident=="G1/S"], rownames(cell_type)[cell_type$ident=="G2/M"])]
regulonAUC_sig_scaled <- pheatmap:::scale_rows(x = regulonAUC_sig)
regulonAUC_sig_scaled <- t(apply(regulonAUC_sig_scaled, 1, function(x) {
  x[x>2.5] <- 2.5
  x[x<(-2.5)] <- (-2.5)
  y <- x
}))

pdf("03-expression/merged/SCENIC_regulonAUC_sig.pdf", width = 6, height = 6)
pheatmap(regulonAUC_sig_scaled, cellwidth = 0.6, show_colnames = F, cluster_rows = F, cluster_cols = F, annotation_col = data.frame(cluster=cell_type$ident, row.names = rownames(cell_type)), gaps_row = 10, gaps_col = 149)
dev.off()

# marker gene information
marker_genes_work <- read.table("03-expression/merged/marker_genes.txt", header = T, sep = "\t", stringsAsFactors = F)
marker_genes_G1_S <- read.table("Data/human/cell_cycle_gene_G1_S.txt", header = F, sep = "\t", stringsAsFactors = F)[,1]
marker_genes_G2_M <- read.table("Data/human/cell_cycle_gene_G2_M.txt", header = F, sep = "\t", stringsAsFactors = F)[,1]
marker_genes_DF <- data.frame(gene=c(marker_genes_work$gene, marker_genes_G1_S, marker_genes_G2_M), 
                              cluster=c(marker_genes_work$cluster, rep("G1/S",length(marker_genes_G1_S)), rep("G2/M", length(marker_genes_G2_M))), 
                              stringsAsFactors = F)
marker_genes_DF$state <- ifelse(marker_genes_DF$gene%in%c(marker_genes_G1_S, marker_genes_G2_M), "list", "novel")
marker_genes_ovlp <- marker_genes_DF$gene[duplicated(marker_genes_DF$gene)]
marker_genes_DF <- marker_genes_DF[ ! duplicated(marker_genes_DF$gene), ]
marker_genes_DF[marker_genes_DF$gene%in%marker_genes_ovlp, "state"] <- "both"
marker_genes_DF$type <- ifelse(marker_genes_DF$gene%in%TF_list, "TF", "gene")
marker_genes_DF <- marker_genes_DF[order(marker_genes_DF$cluster, marker_genes_DF$state), ]
#write.table(file = "03-expression/merged/SCENIC_marker_genes.txt", x = marker_genes_DF, row.names = F, col.names = T, quote = F, sep = "\t")

# detail for specific modules
do_module_table <- function(name, search_mode="TF", add_intra=T, add_inter=T) {
  # base
  ###md <- melt(regulons[name])[,2:1] # regulons was added self TF, which may not in the module
  if(search_mode=="TF") {
    md <- regulonTargetsInfo[regulonTargetsInfo$TF%in%name, 1:2]
  } else if(search_mode=="target") {
    md <- regulonTargetsInfo[regulonTargetsInfo$gene%in%name, 1:2]
  } else if(search_mode=="regulation") {
    md <- regulonTargetsInfo[regulonTargetsInfo$TF%in%name & regulonTargetsInfo$gene%in%name, 1:2]
  } else if(search_mode=="all") {
    md <- regulonTargetsInfo[, 1:2]
  }
  colnames(md) <- c("TF", "target")
  md$group <- "base"
  
  if(add_intra) {
  # add intra-regulon regulation
  md_genes <- unique(c(md$TF, md$target))
  md_intra_added <- regulonTargetsInfo[regulonTargetsInfo$TF%in%md_genes & regulonTargetsInfo$gene%in%md_genes, 1:2]
  colnames(md_intra_added) <- c("TF", "target")
  md_intra_added$group <- "intra"
  md_intra_added <- md_intra_added[md_intra_added$TF!=unique(md$TF), ]
  md <- rbind(md, md_intra_added)
  }
  
  if(add_inter) {
  # add inter-regulon regulation
  md_inter_added <- regulonTargetsInfo[regulonTargetsInfo$TF%in%md_genes | regulonTargetsInfo$gene%in%md_genes, 1:2]
  colnames(md_inter_added) <- c("TF", "target")
  md_inter_added <- setdiff(md_inter_added, md[, 1:2])
  md_inter_added$group <- "inter"
  md <- rbind(md, md_inter_added)
  }

  # more annotation
  md$regulator <- "TF"
  md$targetInTF <- md$target%in%allTFs
  md$interaction <- "pp"
  md$spearman_cor <- apply(md, 1, function(x) { y <- corrMat[x[1], x[2]]; return(y) })
  md$spearman_p <- apply(md, 1, function(x) { y <- cor.test(exprMatr[x[1],], exprMatr[x[2],], method = "spearman"); z <- y$p.value; return(z) })
  md$effect <- apply(md[,7:8], 1, function(x) { if(x[1]>0.03 & x[2]<0.05) { y <- "positive" } else if(x[1]<(-0.03) & x[2]<0.05) { y <- "negative" } else{ y <- "na" }; return(y) })
  md <- md[,c(1,6,2:5,7:9)]
  return(md)
}

# module 1
rg_id <- rownames(head(subset(regulonAUC_csp, ratio<1),1))
rg_tf <- gsub(pattern = " \\(.*", replacement = "", rg_id)
md <- do_module_table(rg_tf)
table(md$TF)
table(md$group)
write.table(file = paste0("03-expression/merged/SCENIC_module_",rg_tf,".txt"), x = md, row.names = F, col.names = T, quote = F, sep = "\t")

# module 2
rg_id <- grep("E2F8",rownames(regulonAUC_csp), value = T)
rg_tf <- gsub(pattern = " \\(.*", replacement = "", rg_id)
md <- do_module_table(rg_tf)
table(md$TF)
table(md$group)
# remove other too big modules
md <- md[- which(md$TF=="E2F2" & (! md$target%in%md$target[md$group=="base"])), ]
write.table(file = paste0("03-expression/merged/SCENIC_module_",rg_tf,".txt"), x = md, row.names = F, col.names = T, quote = F, sep = "\t")

# module 3
rg_id <- grep("E2F3",rownames(regulonAUC_csp), value = T)
rg_tf <- gsub(pattern = " \\(.*", replacement = "", rg_id)
md <- do_module_table(rg_tf)
table(md$TF)
table(md$group)
# remove other too big modules
md <- md[- which(md$TF=="E2F2" & (! md$target%in%md$target[md$group=="base"])), ]
write.table(file = paste0("03-expression/merged/SCENIC_module_",rg_tf,".txt"), x = md, row.names = F, col.names = T, quote = F, sep = "\t")

# module all
rg_tf <- "all"
md <- do_module_table("all", search_mode = "all", add_intra = F, add_inter = F)
dim(md)
table(md$TF)
table(md$group)
write.table(file = paste0("03-expression/merged/SCENIC_module_",rg_tf,".txt"), x = md, row.names = F, col.names = T, quote = F, sep = "\t")

# Spearman correlation issues
g1 <- "E2F7" #"E2F7"
g2 <- "E2F1" #"E2F3"
test_DF <- data.frame(gene1=exprMatr_logCPM[g1, ], gene2=exprMatr_logCPM[g2, ], cell_type=cell_type[colnames(exprMatr),"ident"], stringsAsFactors = F)
tmp <- cor.test(test_DF$gene1, test_DF$gene2, method = "pearson")
data.frame(cor=tmp$estimate, p=tmp$p.value, sig=tmp$p.value<0.05, row.names = "Pearson")
tmp <- cor.test(test_DF$gene1, test_DF$gene2, method = "spearman")
data.frame(cor=tmp$estimate, p=tmp$p.value, sig=tmp$p.value<0.05, row.names = "Spearman")

tmp <- cor.test(test_DF[test_DF$cell_type=="G1/S", ]$gene1, test_DF[test_DF$cell_type=="G1/S", ]$gene2, method = "pearson")
data.frame(cor=tmp$estimate, p=tmp$p.value, sig=tmp$p.value<0.05, row.names = "Pearson")
tmp <- cor.test(test_DF[test_DF$cell_type=="G2/M", ]$gene1, test_DF[test_DF$cell_type=="G2/M", ]$gene2, method = "spearman")
data.frame(cor=tmp$estimate, p=tmp$p.value, sig=tmp$p.value<0.05, row.names = "Spearman")

pdf("03-expression/merged/SCENIC_cor_issue.pdf", width = 6, height = 6, useDingbats = F)
# full
ggplot(test_DF, aes(x=gene1, y=gene2, col=cell_type)) + geom_point(alpha=0.3) + geom_jitter() +
  ggtitle("UMI count") + xlim(c(0,max(test_DF$gene1))) + ylim(c(0,max(test_DF$gene2))) + xlab(g1) + ylab(g2)
# split
ggplot(test_DF[test_DF$cell_type=="G1/S", ], aes(x=gene1, y=gene2, col=cell_type)) + geom_point(alpha=0.3) + geom_jitter() +
  ggtitle("UMI count") + xlim(c(0,max(test_DF$gene1))) + ylim(c(0,max(test_DF$gene2))) + xlab(g1) + ylab(g2) + scale_color_manual(values = scales::hue_pal()(2)[1])
ggplot(test_DF[test_DF$cell_type=="G2/M", ], aes(x=gene1, y=gene2, col=cell_type)) + geom_point(alpha=0.3) + geom_jitter() +
  ggtitle("UMI count") + xlim(c(0,max(test_DF$gene1))) + ylim(c(0,max(test_DF$gene2))) + xlab(g1) + ylab(g2) + scale_color_manual(values = scales::hue_pal()(2)[2])
dev.off()

# mask 0 values
ggplot(test_DF[test_DF$gene1 != 0 & test_DF$gene2 != 0, ], aes(x=gene1, y=gene2, col=cell_type)) + geom_point(alpha=0.3) + geom_jitter() +
  ggtitle("UMI count") + xlim(c(0,max(test_DF$gene1))) + ylim(c(0,max(test_DF$gene2))) + xlab(g1) + ylab(g2)
# split
ggplot(test_DF[test_DF$cell_type=="G1/S" & test_DF$gene1 != 0 & test_DF$gene2 != 0, ], aes(x=gene1, y=gene2, col=cell_type)) + geom_point(alpha=0.3) + geom_jitter() +
  ggtitle("UMI count") + xlim(c(0,max(test_DF$gene1))) + ylim(c(0,max(test_DF$gene2))) + xlab(g1) + ylab(g2) + scale_color_manual(values = scales::hue_pal()(2)[1])
ggplot(test_DF[test_DF$cell_type=="G2/M" & test_DF$gene1 != 0 & test_DF$gene2 != 0, ], aes(x=gene1, y=gene2, col=cell_type)) + geom_point(alpha=0.3) + geom_jitter() +
  ggtitle("UMI count") + xlim(c(0,max(test_DF$gene1))) + ylim(c(0,max(test_DF$gene2))) + xlab(g1) + ylab(g2) + scale_color_manual(values = scales::hue_pal()(2)[2])


# Evidence for E2F7
evd_E2F7 <- read.table(file = "04-GRN/human/00-evidence/01_Westendorp_NAR_2012/network.txt", header = F, sep = "\t", stringsAsFactors = F)
unlist(apply(evd_E2F7, 1, function(x) {
  y <- corrMat[x[1],x[2]]
  return(y)
}))

# the relationship between marker genes
rg_tf <- "marker"
geneList_can <- marker_genes_DF$gene
geneList_DF <- subset(regulonTargetsInfo, gene%in%geneList_can)[, 1:2]  # XXX
geneList_DF$group1 <- marker_genes_DF[match(geneList_DF$TF,marker_genes_DF$gene), "cluster"]
geneList_DF$group2 <- marker_genes_DF[match(geneList_DF$gene,marker_genes_DF$gene), "cluster"]
geneList_DF$regulator <- "TF"
write.table(file = paste0("03-expression/merged/SCENIC_module_",rg_tf,".txt"), x = geneList_DF, row.names = F, col.names = T, quote = F, sep = "\t")

# master regulator for marker genes
# G1/S
geneList_can <- marker_genes_work[marker_genes_work$cluster=="G1/S", "gene"]
test_M <- length(geneList_can)
geneList_DF <- subset(regulonTargetsInfo, gene%in%geneList_can)[, 1:2]  # XXX
test_m <- as.numeric(table(geneList_DF$TF))
test_n <- as.numeric(table(regulonTargetsInfo$TF)[names(table(geneList_DF$TF))])
test_N <- length(unique(regulonTargetsInfo$gene))
test_DF <- data.frame(N=test_N, n=test_n, M=test_M, m=test_m, row.names = names(table(geneList_DF$TF)))

test_DF$p <- apply(test_DF, 1, function(x) {
  y <- phyper(x[4], x[2], x[1]-x[2], x[3], lower.tail = F)
})
test_DF <- test_DF[order(test_DF$p),]
test_DF$q <- p.adjust(p = test_DF$p, method = "BH")
test_DF[test_DF$m>3 & test_DF$q<0.05,]

# the relationship between two modules
load(file="03-expression/merged/GENIE3_corrMat.RData")
dim(corrMat)

{
  geneList_can <- c("E2F3", "KDM7A")
  geneList_DF <- subset(regulonTargetsInfo, TF%in%geneList_can & directAnnot=="TRUE")[, 1:2]  # XXX
  geneList_DF$group <- geneList_DF$TF
  # one gene may be regulated by multiple TFs
  geneList_dup <- unlist(geneList_DF[duplicated(geneList_DF$gene), "gene"])
  geneList_DF <- geneList_DF[! geneList_DF$gene%in%geneList_dup, ]
  geneList_DF <- rbind(geneList_DF, data.frame(TF=rep("mixed",length(geneList_dup)), gene=geneList_dup, group=rep("mixed",length(geneList_dup)), stringsAsFactors = F))
  # add the TF themselves to the genes
  geneList_TF_added <- setdiff(geneList_can, geneList_DF$gene)
  geneList_DF <- rbind(geneList_DF, data.frame(TF=geneList_TF_added, gene=geneList_TF_added, group=geneList_TF_added, stringsAsFactors = F))
  
  geneList_DF <- as.data.frame(geneList_DF)
  geneList_DF <- geneList_DF[order(geneList_DF$group), ]
  rownames(geneList_DF) <- geneList_DF$gene
  geneList <- geneList_DF$gene
  length(geneList)
  corrMat_sub <- corrMat[geneList,geneList]
  dim(corrMat_sub)
  
  ds <- as.dist(1-corrMat_sub)
  hc <- hclust(ds, method = "ward.D2")
  corrMat_sub[corrMat_sub>0.3] <- 0.3
  corrMat_sub[corrMat_sub<(-0.3)] <- (-0.3)
  #corrMat_sub[abs(corrMat_sub)<0.3] <- 0
}

pdf("03-expression/merged/SCENIC_correlation.pdf", width = 6, height = 6, useDingbats = F, onefile = F)
library("pheatmap")
pheatmap(mat = corrMat_sub, cluster_rows = hc, cluster_cols = hc, annotation_row = geneList_DF[, c(-1,-2), drop=F], annotation_col = geneList_DF[, c(-1,-2), drop=F], show_rownames = F, show_colnames = F, color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# find the highly correlated gene set 
gset1 <- names(which(cutree(tree = hc, h = 2.2)==5))
length(gset1)
write.table(x = gset1, file = "03-expression/merged/SCENIC_gset1.txt", row.names = F, col.names = F, quote = F, sep = "\t")
gset2 <- names(which(cutree(tree = hc, h = 2.2)==6))
length(gset2)
write.table(x = gset2, file = "03-expression/merged/SCENIC_gset2.txt", row.names = F, col.names = F, quote = F, sep = "\t")

gset_all <- unique(c("E2F3","KDM7A",gset1,gset2))
md <- do_module_table(name = gset_all, search_mode = "regulation", add_intra = F, add_inter = F)

# gene module detection
nt_genes <- unique(regulonTargetsInfo$gene)
corrMat_sub <- corrMat[nt_genes, nt_genes]
dim(corrMat_sub)
ds <- as.dist(1-corrMat_sub)
hc <- hclust(ds, method = "average")
geneList_DF <- as.data.frame(regulonTargetsInfo[, 2:1])
geneList_dup <- geneList_DF[duplicated(geneList_DF$gene), 1]
geneList_DF <- geneList_DF[ ! duplicated(geneList_DF$gene), ]
geneList_DF[geneList_DF$gene%in%geneList_dup, 2] <- "multiple"
rownames(geneList_DF) <- geneList_DF[,1]
geneList_DF <- geneList_DF[,-1, drop=F]
colnames(geneList_DF) <- "group"

#pdf("03-expression/merged/SCENIC_correlation_sub.pdf", width = 6, height = 6, useDingbats = F)
#tmp_MT <- matrix(1:nrow(corrMat_sub),nrow = 1)
#colnames(tmp_MT) <- nt_genes
#pheatmap(mat = tmp_MT, cluster_rows = F, cluster_cols = hc, annotation_col = geneList_DF, show_rownames = F, show_colnames = F, color = colorRampPalette(c("blue", "white", "red"))(100), cellwidth = 0.2,annotation_legend = F)
#dev.off()








# PCA-based t-SNE
set.seed(123)
tsneAUC <- Rtsne::Rtsne(t(getAuc(regulonAUC_subset)), initial_dims=10, perplexity=10)
rownames(tsneAUC$Y) <- colnames(regulonAUC_subset)
colnames(tsneAUC$Y) <- c("t-SNE1", "t-SNE2")
save(tsneAUC, file="03-expression/merged/SCENIC_tsneRegulonAUC_PCA.RData")

# Alternative: Distance-based t-SNE:
corDist <- as.dist(1-cor(getAuc(regulonAUC_subset)))
set.seed(123)
tsneAUC <- Rtsne::Rtsne(corDist, is_distance=TRUE, perplexity=10)
rownames(tsneAUC$Y) <- labels(corDist)
colnames(tsneAUC$Y) <- c("t-SNE1", "t-SNE2")
save(tsneAUC, file="03-expression/merged/SCENIC_tsneRegulonAUC_Dist.RData")

tSNE <- tsneAUC$Y
par(mfrow=c(1,2))

# Number of genes detected:
nGenesPerCell <- apply(exprMatr, 2, function(x) sum(x>0))
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=F,include.lowest=T))], names(nGenesPerCell))

plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")

# Other known properties:
#for(varName in names(colVars))
#{
#  cellColor <- setNames(colVars[[varName]][cellInfo[,varName]], rownames(cellInfo))
#  plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
#}

#par(mfrow=c(1,4))

# tSNE (colored by number of genes detected per cell)
#plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")
#plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
#plot.new(); plot.new()

# Plot module activity, thresholds & assignment:
#cells_AUCellThresholds <- plot_aucTsne(tSNE=tSNE, exprMat=exprMatr, regulonAUC=regulonAUC, alphaOff=0.1)

