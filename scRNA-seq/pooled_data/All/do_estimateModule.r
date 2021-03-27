
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("parallel")
library("ggplot2")
suppressMessages(library("cowplot"))

samplingPos <- "."
OUT <- paste0("03-expression/merged/estimateModule/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/estimateModule.RData"))

geneSet_LS <- lapply(c("H", paste0("C", 1:7)), function(id) {
  y <- msigdbr::msigdbr(species = "Homo sapiens", category = id)[, c("gs_name", "gene_symbol")]
  y$category <- id
  return(y)
})
geneSet_DF <- do.call("rbind", geneSet_LS)
geneSet_DF$gs_cname <- paste(geneSet_DF$category, geneSet_DF$gs_name, sep = ".")

# gene module
module_genes <- read.table(file = "03-expression/merged/geneModule/geneModule_genes.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(module_genes) <- c("mdid", "gene")
module_genes$mdid <- factor(module_genes$mdid, levels = unique(module_genes$mdid))
module_genes$group <- "GeACT"

module_size <- as.data.frame(table(module_genes$mdid))
colnames(module_size) <- c("mdid", "size")

# random module
random_genes_LS <- lapply(1:1000, function(x) {
  set.seed(x)
  y <- module_genes
  y$gene <- sample(module_genes$gene, nrow(module_genes))
  y$group <- paste0("Random-", x)
  return(y)
})

# combine
com_genes_LS <- c(random_genes_LS, list(module_genes))
names(com_genes_LS) <- c(paste0("Random-", 1:1000), "GeACT")

cl <- makeCluster(50, type = "FORK")
comStat_LS <- parLapply(cl, names(com_genes_LS), function(i) {
  cat(">", i, "\n")
  geneToSet <- merge(x = com_genes_LS[[i]], y = geneSet_DF, by.x = "gene", by.y = "gene_symbol", sort = F)
  geneToSet_LS <- split(geneToSet, geneToSet$mdid)
  
  inGeneSetStat_LS <- lapply(names(geneToSet_LS), function(x) {
    inGeneSet_stat <- table(geneToSet_LS[[x]]$gs_cname)
    inGeneSet_max_num <- max(inGeneSet_stat)
    inGeneSet_max_idx <- which(inGeneSet_stat == inGeneSet_max_num)
    inGeneSet_max_name <- ifelse(length(inGeneSet_max_idx) <= 3, inGeneSet_max_id <- paste(names(inGeneSet_stat)[inGeneSet_max_idx], collapse = ", "), "Many")
    y <- data.frame(mdid = x, inGeneSet_max_num, inGeneSet_max_name, stringsAsFactors = F)
    return(y)
  })
  inGeneSetStat_DF <- do.call("rbind", inGeneSetStat_LS)
  y <- merge(inGeneSetStat_DF, module_size, by = "mdid", sort = F)
  y$inGeneSet_max_ratio <- y$inGeneSet_max_num / y$size
  y$group <- i
  return(y)
})
stopCluster(cl); rm(cl)
comStat_DF <- do.call("rbind", comStat_LS)
comStat_DF$type <- gsub("-.*", "", comStat_DF$group)
comStat_DF$type <- factor(comStat_DF$type, levels = c("Random", "GeACT"))

#ggplot(comStat_DF, aes(x = type, y = inGeneSet_max_ratio, fill = type, color = type)) + geom_violin(show.legend = F)

# average (random)
random_ratio_MT <- do.call("cbind", lapply(comStat_LS[1:1000], "[", "inGeneSet_max_ratio"))
random_ratio_DF <- data.frame(type = "Random", inGeneSet_max_ratio = rowMeans(random_ratio_MT), stringsAsFactors = F)
#
ratio_DF <- rbind(random_ratio_DF, subset(comStat_DF, type == "GeACT", c("type", "inGeneSet_max_ratio")))
ratio_DF$type <- factor(ratio_DF$type, levels = c("Random", "GeACT"))
wilcox.test(subset(ratio_DF, type == "Random", "inGeneSet_max_ratio", drop = T), subset(ratio_DF, type == "GeACT", "inGeneSet_max_ratio", drop = T), alternative = "less")

pdf(file = file.path(OUT, "geneInGeneSet_max_ratio.pdf"), width = 4, height = 4)

ymax <- 1.1
sigBar <- data.frame(type = rep(levels(ratio_DF$type), each = 2), inGeneSet_max_ratio = ymax * c(0.9, 0.925, 0.925, 0.9), stringsAsFactors = F)
ggplot(ratio_DF, aes(x = type, y = inGeneSet_max_ratio)) + geom_violin(aes(fill = type, color = type), show.legend = F) + 
  stat_summary(geom = "point", fun = "median", colour = "white", size = 3, show.legend = F) + 
  geom_path(data = sigBar, group = 1) + annotate(geom = "text", x = 1.5, y = ymax * (1 - 0.05), label = "***", col="red", size = 5.5) + 
  scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 1, by = 0.2)) + 
  scale_fill_manual(values = c("grey75", "cornflowerblue")) + 
  scale_color_manual(values = c("grey75", "cornflowerblue")) + 
  xlab(NULL) + ylab("Ratio of genes in the same gene set") + 
  theme(aspect.ratio = 1)

dev.off()

save.image(file = paste0(OUT, "/estimateModule.RData"))
