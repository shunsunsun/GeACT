# harmony
suppressMessages(library("arrow"))
suppressMessages(library("Seurat"))
suppressMessages(library("harmony"))
source("../../scripts/cluster_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/cellAlign/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/cellAlign.RData"))

# 1. preprocess ----
# load expression matrix
expr_data <- read_feather("03-expression/merged/filtering/UMIcount_filtered.feather")
expr_data <- as.data.frame(expr_data)
expr_data_gene <- read.table("03-expression/merged/filtering/UMIcount_filtered.gene", header = F, sep = "\t", stringsAsFactors = F)
rownames(expr_data) <- expr_data_gene$V1

# load cell meta
cellMeta <- read.table("cell_metatable_filtered_plus.txt", header = T, sep = "\t", stringsAsFactors = F)
all(colnames(expr_data) == cellMeta$cell)

# 2. harmony ----
# norm
expr_norm <- sweep(expr_data, 2, colSums(expr_data), "/") * 1e6
expr_norm_MT <- as.matrix(expr_norm)

# alignment
harmony_embeddings <- HarmonyMatrix(expr_norm_MT, cellMeta, "stage")
rownames(harmony_embeddings) <- cellMeta$cell
colnames(harmony_embeddings) <- paste0("PC", 1:ncol(harmony_embeddings))

# 3. Seurat ----
# Initialize the Seurat object with the raw (non-normalized data).
expr <- CreateSeuratObject(raw.data = expr_data, min.cells = 3, min.genes = 500, project = samplingPos, names.delim = "/")
dim(expr@raw.data)

# add meta
expr@meta.data <- cbind(expr@meta.data, cellMeta[match(rownames(expr@meta.data), cellMeta$cell), c("stage", "ident")])
expr@meta.data$cluster <- expr@meta.data$ident
mito.genes <- grep(pattern = "^MT-", x = rownames(expr@raw.data))
length(mito.genes)
percent.mito <- Matrix::colSums(expr@raw.data[mito.genes, ])/Matrix::colSums(expr@raw.data)
expr <- do_addMeta(expr)

expr <- FilterCells(object = expr, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, Inf))
expr <- NormalizeData(object = expr, normalization.method = "LogNormalize", scale.factor = 10000)
expr <- FindVariableGenes(object = expr, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.25, x.high.cutoff = 5, y.cutoff = 0.5, do.plot = F)
expr <- ScaleData(object = expr, vars.to.regress = c("nUMI", "percent.mito"), num.cores = 20, do.par = T)

# run PCA
expr <- RunPCA(object = expr, pc.genes = expr@var.genes, pcs.compute = 20, do.print = F)

# create harmony embeddings
expr@dr$harmony <- expr@dr$pca
expr@dr$harmony@cell.embeddings <- harmony_embeddings
expr@dr$harmony@gene.loadings <- matrix()
expr@dr$harmony@sdev <- numeric()

# run tSNE/UMAP (based on Harmony)
expr <- RunTSNE(object = expr, reduction.use = "harmony", dims.use = 1:20, nthreads = 20, do.fast = T)
expr <- RunUMAP(object = expr, reduction.use = "harmony", dims.use = 1:20, min_dist = 1)

# plot
# DimPlot(object = expr, reduction.use = "pca", pt.size = 2, do.label = T, no.legend = T, plot.title = "PCA", group.by = "stage")
# DimPlot(object = expr, reduction.use = "pca", pt.size = 2, do.label = T, no.legend = T, plot.title = "PCA", group.by = "cluster")
# DimPlot(object = expr, reduction.use = "harmony", pt.size = 2, do.label = T, no.legend = T, plot.title = "Harmony", group.by = "stage")
# DimPlot(object = expr, reduction.use = "harmony", pt.size = 2, do.label = T, no.legend = T, plot.title = "Harmony", group.by = "cluster")
# DimPlot(object = expr, reduction.use = "tsne", pt.size = 2, do.label = T, no.legend = T, plot.title = "tSNE", group.by = "stage")
# DimPlot(object = expr, reduction.use = "tsne", pt.size = 2, do.label = T, no.legend = T, plot.title = "tSNE", group.by = "cluster")
# DimPlot(object = expr, reduction.use = "umap", pt.size = 2, do.label = T, no.legend = T, plot.title = "tSNE", group.by = "stage")
# DimPlot(object = expr, reduction.use = "umap", pt.size = 2, do.label = T, no.legend = T, plot.title = "tSNE", group.by = "cluster")

cellMetaData <- merge(expr@meta.data, cbind(expr@dr$tsne@cell.embeddings, expr@dr$umap@cell.embeddings), by = 0, sort = F)
colnames(cellMetaData)[1] <- "cell"
colnames(cellMetaData)[grep("^UMAP", colnames(cellMetaData))] <- c("UMAP_1", "UMAP_2")

# 4. save and write meta table ----
write.table(x = cellMetaData, file = paste0(OUT, "/Seurat_metaData.txt"), row.names = F, col.names = T, quote = F,sep = "\t")

save.image(file = paste0(OUT, "/cellAlign.RData"))
