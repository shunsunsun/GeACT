cc.preprocess <- function(expr_data, cellStat, project = "TestProject", min.cells = 3, min.genes = 500) {
  # Initialize the Seurat object with the raw (non-normalized data).
  expr <- CreateSeuratObject(raw.data = expr_data, min.cells = min.cells, min.genes = min.genes, project = project, names.delim = "/")

  # QC and selecting cells for further analysis
  # calculate the percent.mito values.
  mito.genes <- grep(pattern = "^MT-", x = rownames(expr@raw.data))
  percent.mito <- Matrix::colSums(expr@raw.data[mito.genes, ])/Matrix::colSums(expr@raw.data)
  expr <- AddMetaData(object = expr, metadata = percent.mito, col.name = "percent.mito")
  
  # add batch info
  if(length(grep("_cell[0-9]+$", colnames(expr@raw.data))) == ncol(expr@raw.data)) {
    batch <- gsub("_cell.*", "", colnames(expr@raw.data))
    names(batch) <- colnames(expr@raw.data)
    bigBatch <- gsub("-[0-9]*_cell.*", "", colnames(expr@raw.data))
    names(bigBatch) <- colnames(expr@raw.data)
    #people <- rep(c("A", "B"), c(792, 2463-792))
    #names(people) <- colnames(expr@raw.data)
    outer_id <- ceiling(as.numeric(gsub(".*cell", "", colnames(expr@raw.data))) / 8); outer_id <- factor(outer_id)
    names(outer_id) <- colnames(expr@raw.data)
    inner_id <- ceiling(as.numeric(gsub(".*cell", "", colnames(expr@raw.data))) %% 8); inner_id[inner_id==0] <- 8; inner_id <- factor(inner_id, levels = 1:8)
    names(inner_id) <- colnames(expr@raw.data)
    bigInner_id <- ifelse(as.numeric(inner_id)<=4, "1-4", "5-8")
    names(bigInner_id) <- colnames(expr@raw.data)
    
    expr <- AddMetaData(object = expr, metadata = batch, col.name = "batch")
    expr <- AddMetaData(object = expr, metadata = bigBatch, col.name = "bigBatch")
    #expr <- AddMetaData(object = expr, metadata = people, col.name = "people")
    expr <- AddMetaData(object = expr, metadata = outer_id, col.name = "outer_id")
    expr <- AddMetaData(object = expr, metadata = inner_id, col.name = "inner_id")
    expr <- AddMetaData(object = expr, metadata = bigInner_id, col.name = "bigInner_id")
  } else {
    warning("The cell ID do not match _cell[0-9], so the batch label will be skipped.")
  }
  
  #VlnPlot(object = expr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.title.use = 16)

  # cell filtering
  expr <- FilterCells(object = expr, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, Inf))

  # Normalizing the data
  expr <- NormalizeData(object = expr, normalization.method = "LogNormalize", scale.factor = 10000)

  # Detection of variable genes across the single cells
  expr <- FindVariableGenes(object = expr, mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 1)
  # length(x = expr@var.genes)

  #Scaling the data and removing unwanted sources of variation
  expr <- ScaleData(object = expr, vars.to.regress = c("nUMI", "percent.mito"), num.cores = 15, do.par = T)
  #expr <- ScaleData(object = expr, vars.to.regress = c("nUMI"), num.cores = 15, do.par = T)

  return(expr)
}

cc.projection <- function(expr, pcs.compute = 100) {
  #Perform linear dimensional reduction
  expr <- RunPCA(object = expr, pc.genes = expr@var.genes, pcs.compute = pcs.compute, do.print = F)

  # Examine and visualize PCA results a few different ways
  #PrintPCA(object = expr, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  #par(oma=c(0,2,0,0))
  #VizPCA(object = expr, pcs.use = 1:2)
  #PCAPlot(object = expr, dim.1 = 1, dim.2 = 2)

  # ProjectPCA scores each gene in the dataset (including genes not included
  # in the PCA) based on their correlation with the calculated components.
  expr <- ProjectPCA(object = expr, do.print = FALSE)
  return(expr)
}

cc.projection.visual <- function(expr, force.recalc = F) {
  # Determine statistically significant principal components
  if(force.recalc | is.null(expr@dr$pca@jackstraw)) {
    expr <- JackStraw(object = expr, num.replicate = 100, num.cores = 15, do.par = T, maxit = 1000)
  }
  JackStrawPlot(object = expr, PCs = 1:20)
  PCElbowPlot(object = expr, num.pc = 30)
  return(expr)
}

cc.clustering <- function(expr, dims_use=1:20, resolution = 0.6) {
  # Cluster the cells
  expr <- FindClusters(object = expr, reduction.type = "pca", dims.use = dims_use,
                       resolution = resolution, print.output = 0, save.SNN = TRUE, temp.file.location = "/tmp/")
  #PrintFindClustersParams(object = expr)
  expr@meta.data$cluster <- expr@meta.data[, grep("res.", colnames(expr@meta.data), fixed = T)]
  #table(expr@meta.data$cluster)

  expr <- RunTSNE(object = expr, dims.use = dims_use, do.fast = TRUE)
  TSNEPlot(object = expr, pt.size = 2, do.label = T, no.legend = T)

  return(expr)
}
