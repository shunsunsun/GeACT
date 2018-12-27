# gene module detection
library("pheatmap")
library("gi")
# source("src/module.r") # >>>

# Load the dataset
projectID <- "geact2"
sessionID <- "5c24e64af20f9c1358adcd59"

expr_data <- t(getExprMatrix(projectID, sessionID))
cellStat <- getMeta(projectID, sessionID)
cellStat$expr.ident <- "Mono"

expr_sub_LS <- md.preprocess(expr_data, cellStat, outdir = "tmp") # cell filtering and split cells by cluster
hc_cor_LS <- md.clustering(expr_sub_LS)  # hierarchy clustering
cluster_table_LS <- md.module(hc_cor_LS, outdir = "tmp") # module detection based on clustering

# save results as RDS
# saveRDS(object = cluster_table_LS, file = "do_module.RDS")
