# gene module detection
library("pheatmap")
library("gi")
source("src/module.r") # >>>

# Load the dataset
projectID <- "geact2"
sessionID <- "5c2445c9f20f9c1358adcb86"

expr_data <- t(getExprMatrix(projectID, sessionID))
cellStat <- getMeta(projectID, sessionID)
cellStat$expr.ident <- "Mono"

expr_sub_LS <- step1_preprocess(expr_data, cellStat, outdir = "tmp") # cell filtering and split cells by cluster
hc_cor_LS <- step2_clustering(expr_sub_LS)  # hierarchy clustering
cluster_table_LS <- step3_module(hc_cor_LS, outdir = "tmp") # module detection based on clustering

saveRDS(object = cluster_table_LS, file = "do_module.RDS")
