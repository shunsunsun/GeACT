# cell classification
library("Seurat")
library("dplyr")
library("Matrix")
library("gi")
# source("src/cluster.r") # >>>

# Load the dataset
projectID <- "geact2"
sessionID <- "5c24e5c9f20f9c1358adcd4d"
expr_data <- t(getExprMatrix(projectID, sessionID))
cellStat <- getMeta(projectID, sessionID)

# running
expr <- cc.preprocess(expr_data, cellStat) # filter genes/cells and normalization
expr <- cc.projection(expr) # run PCA and show significant dimensions
expr <- cc.clustering(expr, dims_use = 1:13, resolution = 0.6) # cluster cells and run tSNE

# save results as RDS
# saveRDS(object = expr, file = "do_cluster.Rds")
