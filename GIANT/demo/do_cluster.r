# cell classification
library(Seurat)
library(dplyr)
library(Matrix)
library("gi")
source("src/cluster.r") # >>>

# Load the dataset
projectID <- "geact2"
sessionID <- "5c2325197e71c30a28f5359b"
expr_data <- t(getExprMatrix(projectID, sessionID))
cellStat <- getMeta(projectID, sessionID)

# running
expr <- step1_preprocess(expr_data, cellStat) # filter genes/cells and normalization
expr <- step2_chooseDim(expr) # run PCA and show significant dimensions
expr <- step3_clustering(expr, dims_use = 1:20, resolution = 0.6) # cluster cells and run tSNE
# save results as RDS
saveRDS(object = expr, file = "do_cluster.Rds")
