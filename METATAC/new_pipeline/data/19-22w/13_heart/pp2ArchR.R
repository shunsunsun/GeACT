#!/usr/bin/env r
suppressMessages({
  library(docopt)
})

# Argument parsing -------------------------------------------------------------
doc <- "Usage: main_ArchR_pipeline.R [-h] [--root ROOT]

-h --help           show help message
--root ROOT         root directory of ATAC analysis [default: /data/Lab/otherwork/GeACT/ATAC]
"

opt <- docopt(doc)

if(opt$help){
  cat(doc)
  q("no")
}

root <- opt$root

suppressMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(gridExtra)
  library(Seurat)
})

# env setup --------------------------------------------------------------------
age <- "19-22w"
tissue <- "13_heart"
.data <- "data"
.lix_script <- "scripts_lix"

organ_wd <- paste(root, .data, age, tissue, sep = '/')
setwd(organ_wd)

lix_script_dir <- paste(root, .lix_script, sep = "/")
source(paste(lix_script_dir, "QC_utils.R", sep = "/"))

# 0 Cell metadata setup --------------------------------------------------------

meta_table <- read.delim("../../../meta/meta_table_heart.txt", row.names = 1,
                         stringsAsFactors = F, check.names = F)
cellMetaData <- read.delim("frag_and_meta/mapStat_human.txt", row.names = 1, check.names = F)
cellMetaData <- cellMetaData[, c(1, 6)]
cellMetaData$plate <- gsub("_[0-9]*-[0-9]$", "", rownames(cellMetaData))
cellMetaData$seqID <- meta_table[cellMetaData$plate, "sequencing_time"]
cellMetaData$individual <- meta_table[cellMetaData$plate, "individual"]
cellMetaData$tissue <- meta_table[cellMetaData$plate, "tissue"]
cellMetaData$samplingPos <- meta_table[cellMetaData$plate, "samplingPos"]
cellMetaData <- cellMetaData[, c(4:7,3,1,2)]
rm(meta_table)

rownames(cellMetaData) = gsub("P[RD]10_HCA_15_","",rownames(cellMetaData))
rownames(cellMetaData) = gsub("P[RD]10_HCA_","",rownames(cellMetaData))
rownames(cellMetaData) = gsub("Plate","",rownames(cellMetaData))

num_frag <- read.delim("frag_and_meta/num_frag_human.txt",header = F,row.names = 1,check.names = F)
num_frag_decon <- read.delim("frag_and_meta/num_frag_decon_human.txt",header = F,row.names = 1,check.names = F)
num_mito_frag_decon <- read.delim("frag_and_meta/num_mito_frag_decon_human.txt",header = F,row.names = 1,check.names = F)
cellMetaData$num_frag <- num_frag[rownames(cellMetaData),1]
cellMetaData$num_frag_decon <- num_frag_decon[rownames(cellMetaData),1]
cellMetaData$num_mito_frag_decon <- num_mito_frag_decon[rownames(cellMetaData),1]
rm(num_frag,num_frag_decon,num_mito_frag_decon)
cellMetaData$num_mito_frag_decon[is.na(cellMetaData$num_mito_frag_decon)] <- 0

cellMetaData$con_rate <- 1 - cellMetaData$num_frag_decon/cellMetaData$num_frag
cellMetaData$mito_ratio <- cellMetaData$num_mito_frag_decon/(cellMetaData$num_frag_decon + cellMetaData$num_mito_frag_decon)



# 2 ArchR process --------------------------------------------------------------
# setup
ArchR_wd <- paste(root, .data, age, tissue, "results/ArchR", sep = '/')
dir.create(ArchR_wd, recursive = T, showWarnings = F)
setwd(ArchR_wd)
addArchRThreads(16)
# addArchRGenome("hg38")

geneAnnotation <- readRDS(paste(root, "database", "annotation", "geneAnnotation.rds", sep = "/"))
genomeAnnotation <- readRDS(paste(root, "database", "annotation", "genomeAnnotation.rds", sep = "/"))
# genomeAnnotation <- getArchRGenome(genomeAnnotation = T)

# 2.1 Creating Arrow Files -----------------------------------------------------
ArrowFiles <- createArrowFiles(
  inputFiles = c("../../frag_and_meta/merge_human_frag_decon.bed.gz"),
  sampleNames = c(paste(age,tissue,sep = '_')),
  minTSS = 0,
  minFrags = 0,
  maxFrags = Inf,
  addTileMat = T,
  addGeneScoreMat = T,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  force = T,
  subThreading = F
)

# Inferring scATAC-seq Doublets with ArchR
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
) 

# 2.2 Creating An ArchRProject -------------------------------------------------
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR_output",
  copyArrows = T,
  geneAnnotation = geneAnnotation, 
  genomeAnnotation = genomeAnnotation, 
  showLogo = F
)

# 2.3 QC and filtering ---------------------------------------------------------
rownames(cellMetaData) <- paste(paste(age, tissue,sep = '_'), rownames(cellMetaData), sep = '#')
cellMetaData <- cbind(cellMetaData, data.frame(proj@cellColData)[rownames(cellMetaData), ])
# this result is no longer equal due to adding chrM in genome annotation
# all.equal(cellMetaData$num_frag_decon, cellMetaData$nFrags)
# cellMetaData <- cellMetaData[, -which(names(cellMetaData) == 'nFrags')]

hist(log10(cellMetaData$Reads),breaks = 20)
hist(cellMetaData$Reads,breaks = 20)
hist(cellMetaData$Aligned_ratio,breaks = 20)
hist(log10(cellMetaData$num_frag_decon),breaks = 20)
hist(cellMetaData$con_rate, breaks = 20)
hist(cellMetaData$mito_ratio, breaks = 20)
hist(cellMetaData$TSSEnrichment,breaks = 20)
hist(cellMetaData$PromoterRatio,breaks = 20)
hist(cellMetaData$NucleosomeRatio,breaks = 20)
hist(cellMetaData$DoubletEnrichment,breaks = 20)

nReads_l = 1e5
nReads_h = 6e6
aligned_l = 0.85
nFrags_l = log10(1e4)
nFrags_h = 5.4
con_rate_h = 0.9
mito_h = 0.1
TSS_l = 5
Promoter_l = 0.1
Doublet_h = 20

keep_DF <- data.frame(
  reads = cellMetaData$Reads >= nReads_l & cellMetaData$Reads <= nReads_h,
  aligned_ratio = cellMetaData$Aligned_ratio >= aligned_l,
  num_frags = log10(cellMetaData$num_frag_decon) >= nFrags_l & log10(cellMetaData$num_frag_decon) <= nFrags_h,
  con_rate = cellMetaData$con_rate <= con_rate_h,
  mito_ratio = cellMetaData$mito_ratio <= mito_h,
  TSSEnrichment = cellMetaData$TSSEnrichment >= TSS_l,
  PromoterRatio = cellMetaData$PromoterRatio >= Promoter_l,
  DoubletEnrichment = cellMetaData$DoubletEnrichment <= Doublet_h
)
keep_DF <- keep_DF*1
rownames(keep_DF) <- rownames(cellMetaData)
sum(rowSums(keep_DF) == ncol(keep_DF))

cellMetaData$keep <- (rowSums(keep_DF)==ncol(keep_DF))

draw_QCplot("../QC_plot.pdf",cellMetaData,keep_DF,nReads_l,nReads_h,aligned_l,nFrags_l,nFrags_h,
            con_rate_h,mito_h,TSS_l,Promoter_l,Doublet_h)

# do the filtering
kept_cell <- rownames(cellMetaData)[cellMetaData$keep]
write.table(cellMetaData,"../cellMetaData.txt",sep = '\t',row.names = T,col.names = T,quote = F)
# write.table(keep_DF,"../keep_DF.txt",sep = '\t',row.names = T,col.names = T,quote = F)
proj <- proj[kept_cell, ]
nCells(proj)
cellMetaData <- cellMetaData[cellMetaData$keep, ]

# add metadata
proj$seqID = cellMetaData$seqID
proj$individual = cellMetaData$individual
proj$tissue = cellMetaData$tissue
proj$samplingPos = cellMetaData$samplingPos
proj$plate = cellMetaData$plate
proj$Reads = cellMetaData$Reads
proj$Aligned_ratio = cellMetaData$Aligned_ratio
proj$num_frag = cellMetaData$num_frag
proj$num_mito_frag_decon = cellMetaData$num_mito_frag_decon
proj$con_rate = cellMetaData$con_rate
proj$mito_ratio = cellMetaData$mito_ratio

# remove cells from FDM
proj <- proj[!is.na(proj$samplingPos), ]

# Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1, p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# 2.4 Saving and Loading an ArchRProject ---------------------------------------
saveArchRProject(ArchRProj = proj, outputDirectory = "ArchR_output", load = F, dropCells = T)

