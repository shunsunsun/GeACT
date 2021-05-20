suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
})

setwd("/data/Lab/ownwork/GeACT/ATAC")

addArchRGenome("hg19")
addArchRThreads(threads = 16)

inputFiles <- getTutorialData("Hematopoiesis")


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = F
)

# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#   LSIMethod = 1
# )

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

proj_s1 <- projHeme1[projHeme1$Sample == "scATAC_BMMC_R1", ]
proj_s1 <- proj_s1[proj_s1$PromoterRatio > 0.1, ]

proj_s1 <- addIterativeLSI(
  ArchRProj = proj_s1,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 1,
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj_s1 <- addGeneScoreMatrix(proj_s1)

getAvailableMatrices(proj_s1)
getAvailableMatrices(projHeme1)



saveArchRProject(proj_s1)

projHeme1 <- loadArchRProject("HemeTutorial/")


a <- 10

test <- function() {
  a <- a + 20
  print(a)
}

test()
a
