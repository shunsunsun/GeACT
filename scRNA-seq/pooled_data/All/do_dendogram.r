# cell classification
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("Seurat")
library("dplyr")
library("Matrix")
library("pheatmap")
library("dendextend")
library("parallel")
source("../../scripts/cellType_tools.r")
#source("../../scripts/cluster_tools.r")

samplingPos <- "."
OUT <- paste0("03-expression/merged/dendogram/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/dendrogram.RData"))

TF <- read.table("../../Data/human_TF.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(TF)
lncRNA <- read.table("../../Data/human_lncRNA.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(lncRNA)
chromtainRM <- read.table("../../Data/human_cr.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(chromtainRM)
cellSurface <- read.table("../../Data/human_cs.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(cellSurface)

clgrp2color <- function(clgrp.input, color.input = c("tomato", "orange", "#1ADBB8", "#00C9C9", "#00BF74", 
                                                     "#8258FA", "#D358F7", "#D0A9F5", "dodgerblue2", "firebrick3")) {
  clgrp.input <- gsub(".*\\.", "", clgrp.input)
  clgrp.input[clgrp.input == "Epithelial"] <- color.input[1]
  clgrp.input[clgrp.input == "Endothelial"] <- color.input[2]
  clgrp.input[clgrp.input == "Smooth muscle"] <- color.input[3]
  clgrp.input[clgrp.input == "Skeletal muscle"] <- color.input[4]
  clgrp.input[clgrp.input == "Fibroblast"] <- color.input[5]
  clgrp.input[clgrp.input == "B"] <- color.input[6]
  clgrp.input[clgrp.input == "T"] <- color.input[7]
  clgrp.input[clgrp.input %in% c("DC", "DC/Macrophage", "Mast", "Neutrophil", "NKT", "Immune")] <- color.input[8]
  clgrp.input[clgrp.input == "CACNA1A"] <- color.input[9]
  clgrp.input[clgrp.input == "Erythrocyte"] <- color.input[10]
  clgrp.input[! clgrp.input %in% color.input] <- "black"
  return(clgrp.input)
}

plot_hc <- function(hnorm.input, dend.input, cex.para = 1, mar.para = NULL, ylim.para, main = NULL) {
  par_old <- par()
  par(cex = cex.para)
  if(! is.null(mar.para)) { par(mar = mar.para) }
  labels_colors(dend.input) <- clgrp2color(hnorm.input$labels[hnorm.input$order])
  plot(dend.input, horiz = T, ylim = ylim.para, main = main)
  par(par_old[c("cex", "mar")])
}

# 1. pre-process ----
# Load gene ID 
geneID <- read.table("~/lustre/06-Human_cell_atlas/Genomes/human/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(geneID)
colnames(geneID) <- c("ensembl", "symbol")

# Load the expr matrix
expr_data <- read.table(file = paste0("03-expression/merged/filtering/", samplingPos, "/UMIcount_filtered.txt"), header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F, comment.char = "")
dim(expr_data)
# norm
expr_data_normed <- sweep(expr_data, 2, colSums(expr_data), "/")

# Load the cell metatable
cellMetaData <- read.table(file = "cell_metatable_filtered.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
dim(cellMetaData)
all(colnames(expr_data) == rownames(cellMetaData))
cellMetaData$ts_ident <- Hmisc::capitalize(paste(cellMetaData$tissue, cellMetaData$ident, sep = "."))
cellMetaData$ts_clgrp <- Hmisc::capitalize(paste(cellMetaData$tissue, ident2clgrp(cellMetaData$ident), sep = "."))
length(unique(cellMetaData$ts_clgrp))
#View(unique(cellMetaData$ts_clgrp))

# avg
expr_data_avg <- sapply(split(rownames(cellMetaData), cellMetaData$ts_clgrp), function(x) { y <- rowMeans(expr_data_normed[, x, drop = F]) })
# class by immune
immune_clgrp <- c("B", "DC", "DC/Macrophage", "Mast", "Neutrophil", "NKT", "T", "Immune")
expr_data_avg_nonImmune <- expr_data_avg[, ! sapply(strsplit(colnames(expr_data_avg), ".", fixed = T), "[", 2) %in% immune_clgrp]
expr_data_avg_Immune <- expr_data_avg[, sapply(strsplit(colnames(expr_data_avg), ".", fixed = T), "[", 2) %in% immune_clgrp]

# 2.1 hc: all clgrp ----

# 1) all genes
dnorm_all <- as.dist(1 - cor(expr_data_avg, method = "pearson"))
hnorm_all <- hclust(dnorm_all, method = "average")
dend_all <- as.dendrogram(hnorm_all)

# 2) TF
dnorm_tf <- as.dist(1 - cor(expr_data_avg[rownames(expr_data_avg) %in% TF$V1, ], method = "pearson"))
hnorm_tf <- hclust(dnorm_tf, method = "average")
dend_tf <- as.dendrogram(hnorm_tf)

# 3) chromatin remodeler
dnorm_cr <- as.dist(1 - cor(expr_data_avg[rownames(expr_data_avg) %in% chromtainRM$V2, ], method = "pearson"))
hnorm_cr <- hclust(dnorm_cr, method = "average")
dend_cr <- as.dendrogram(hnorm_cr)

# 4) lncRNA
dnorm_lnc <- as.dist(1 - cor(expr_data_avg[rownames(expr_data_avg) %in% lncRNA$V1, ], method = "pearson"))
hnorm_lnc <- hclust(dnorm_lnc, method = "average")
dend_lnc <- as.dendrogram(hnorm_lnc)

# 5) cell surface marker
dnorm_cs <- as.dist(1 - cor(expr_data_avg[rownames(expr_data_avg) %in% chromtainRM$V2, ], method = "pearson"))
hnorm_cs <- hclust(dnorm_cs, method = "average")
dend_cs <- as.dendrogram(hnorm_cs)

# N) module map
# module_map <- read.table(paste0("03-expression/merged/geneModule/", samplingPos, "/module_map.txt"), header = T, sep = "\t", row.names = 1, check.names = F)
# module_map_norm <- sweep(module_map, 1, rowSums(module_map, na.rm = T), "/")
# dnorm_mdm <- dist(module_map_norm)
# hnorm_mdm <- hclust(dnorm_mdm, method = "ward.D2", members = NULL)
# plot(as.dendrogram(hnorm_mdm), cex=.5, horiz=T)

pdf(paste0(OUT, "/hc_all.pdf"), width = 3, height = 18)

plot_hc(hnorm_all, dend_all, cex.para = 0.8, mar.para = c(2,1,0,14) + 0.1, ylim.para = c(3, 130))
plot_hc(hnorm_tf, dend_tf, cex.para = 0.8, mar.para = c(2,1,0,14) + 0.1, ylim.para = c(3, 130))
plot_hc(hnorm_cr, dend_cr, cex.para = 0.8, mar.para = c(2,1,0,14) + 0.1, ylim.para = c(3, 130))
plot_hc(hnorm_lnc, dend_lnc, cex.para = 0.8, mar.para = c(2,1,0,14) + 0.1, ylim.para = c(3, 130))
plot_hc(hnorm_cs, dend_cs, cex.para = 0.8, mar.para = c(2,1,0,14) + 0.1, ylim.para = c(3, 130))

dev.off()

do_xxx <- function() {

# compare
dend1 <- as.dendrogram(hnorm_all)
dend2 <- as.dendrogram(hnorm_tf)
dend3 <- as.dendrogram(hnorm_cr)
dend4 <- as.dendrogram(hnorm_lnc)
#dend5 <- as.dendrogram(hnorm_mdm)

dend12 <- dendlist(dend1, dend2)
dend13 <- dendlist(dend1, dend3)
dend14 <- dendlist(dend1, dend4)
#dend15 <- dendlist(dend1, dend5)

enta12 <- entanglement(untangle(dend12, method = "step2side"))
#enta13 <- entanglement(untangle(dend13, method = "step2side"))
#enta14 <- entanglement(untangle(dend14, method = "step2side"))
#enta15 <- entanglement(untangle(dend15, method = "step2side"))

# random (TF number)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = c("expr_data_all_grp", "dendlist", "dend1", "%>%", "untangle", "entanglement"))
# do_rnd_entanglement <- function(genenum = 1000, repnum = 1000) {
#   clusterExport(cl, varlist = c("genenum"), envir = environment())
#   rnd_value <- parSapply(cl, 1:repnum, function(i) {
#     print(i)
#     set.seed(i)
#     expr_data_rnd_grp <- expr_data_all_grp[, sample(ncol(expr_data_all_grp), genenum)]
#     expr_data_rnd_grp_norm <- sweep(expr_data_rnd_grp, 1, rowSums(expr_data_rnd_grp), "/")
#     dnorm_rnd <- dist(expr_data_rnd_grp_norm)
#     hnorm_rnd <- hclust(dnorm_rnd, method = "ward.D2", members = NULL)
#     
#     dends <- dendlist(dend1, as.dendrogram(hnorm_rnd))
#     y <- dends %>% untangle(method = "step2side")
#     out <- entanglement(y)
#     return(out)
#   })
#   return(rnd_value)
# }
# rnd_LS <- list()
# for(num in c(50, 100, nrow(expr_data_tf), nrow(expr_data_cr), nrow(expr_data_lnc), nrow(expr_data))) {
#   print(num)
#   rnd_LS[[as.character(num)]] <- data.frame(gene_num = num, ent = do_rnd_entanglement(genenum = num, repnum = 1000))
# }
# stopCluster(cl); rm(cl)
# rnd_DF <- do.call("rbind", rnd_LS)

# plot
pdf(paste0(OUT, "/dendrogram.pdf"), width = 6, height = 6)

# ggplot(rnd_DF, aes(x = gene_num, y = ent)) + 
#   stat_summary(fun.y=mean,geom="line", show.legend = F) + 
#   stat_summary(fun.data=mean_se, geom = "errorbar", show.legend = F) + 
#   geom_point(data = data.frame(gene_num = nrow(expr_data_tf), ent = enta12), color = scales::hue_pal()(3)[1]) + 
#   geom_point(data = data.frame(gene_num = nrow(expr_data_cr), ent = enta13), color = scales::hue_pal()(3)[2]) + 
#   geom_point(data = data.frame(gene_num = nrow(expr_data_lnc), ent = enta14), color = scales::hue_pal()(3)[3]) + 
#   geom_text(data = data.frame(gene_num = nrow(expr_data_tf), ent = enta12), label = "TFs", vjust = -1) + 
#   geom_text(data = data.frame(gene_num = nrow(expr_data_cr), ent = enta13), label = "Chromtain remodelers", vjust = 2) + 
#   geom_text(data = data.frame(gene_num = nrow(expr_data_lnc), ent = enta14), label = "LncRNAs", vjust = 2) + 
#   xlab("Gene number") + ylab("Entanglement") + scale_x_log10()

#dend_diff(dend1, dend2)
dl12 <- dend12 %>% untangle(method = "step2side")
tanglegram(dl12, lab.cex = 1.2, cex_main = 1, cex_main_left = 1.25, cex_main_right = 1.25, 
           lwd = 1, edge.lwd = 1, common_subtrees_color_branches = TRUE,common_subtrees_color_lines = TRUE,
           columns_width = c(15, 10, 15), main_left = paste0("All genes\n", nrow(expr_data)), main_right = paste0("Transcription factors\n", sum(rownames(expr_data_avg) %in% TF$V1)),
           margin_inner= 6.5, main = paste("Entanglement =", round(entanglement(dl12), 3)))

# dl13 <- dend13 %>% untangle(method = "step2side")
# tanglegram(dl13, lab.cex = 1.2, cex_main = 1, cex_main_left = 1.25, cex_main_right = 1.25, 
#            lwd = 1, edge.lwd = 1, common_subtrees_color_branches = TRUE,common_subtrees_color_lines = TRUE,
#            columns_width = c(15, 10, 15), main_left = paste0("All genes\n", nrow(expr_data)), main_right = paste0("Chromtain remodelers\n", nrow(expr_data_cr)),
#            margin_inner= 6.5, main = paste("Entanglement =", round(entanglement(dl13), 3)))
# 
# dl14 <- dend14 %>% untangle(method = "step2side")
# tanglegram(dl14, lab.cex = 1.2, cex_main = 1, cex_main_left = 1.25, cex_main_right = 1.25, 
#            lwd = 1, edge.lwd = 1, common_subtrees_color_branches = TRUE,common_subtrees_color_lines = TRUE,
#            columns_width = c(15, 10, 15), main_left = paste0("All genes\n", nrow(expr_data)), main_right = paste0("LncRNAs\n", nrow(expr_data_lnc)), 
#            margin_inner= 6.5, main = paste("Entanglement =", round(entanglement(dl14), 3)))
# 
# dl15 <- dend15 %>% untangle(method = "step2side")
# tanglegram(dl15, lab.cex = 1.2, cex_main = 1, cex_main_left = 1.25, cex_main_right = 1.25, 
#            lwd = 1, edge.lwd = 1, common_subtrees_color_branches = TRUE,common_subtrees_color_lines = TRUE,
#            columns_width = c(15, 10, 15), main_left = paste0("All genes\n", nrow(expr_data)), main_right = paste0("Modules\n", nrow(t(module_map))), 
#            margin_inner= 6.5, main = paste("Entanglement =", round(entanglement(dl15), 3)))

# tf cell type specific
plot.new()
p <- pheatmap(tmp_t_sorted, scale = "none", color = colorRampPalette(c("grey", "navy"))(100), angle_col = 90, 
         show_rownames = F, cluster_rows = F, cluster_cols = F, gaps_row = tmp_n, silent = T)
print(p)

plot.new()
p <- pheatmap(gene_cluster_binned_sorted, scale = "none", color = colorRampPalette(c("black", "limegreen"))(100), legend = T, angle_col = 90, 
              show_rownames = F, cluster_rows = F, cluster_cols = F, gaps_row = sum(rowSums(gene_cluster_binned>0) >= 10), silent = T)
print(p)

dev.off()

}

# 2.2 hc: nonImmune clgrp ----

# 1) all genes
dnorm_all <- as.dist(1 - cor(expr_data_avg_nonImmune, method = "pearson"))
hnorm_all <- hclust(dnorm_all, method = "average")
dend_all <- as.dendrogram(hnorm_all)

# 2) TF
dnorm_tf <- as.dist(1 - cor(expr_data_avg_nonImmune[rownames(expr_data_avg_nonImmune) %in% TF$V1, ], method = "pearson"))
hnorm_tf <- hclust(dnorm_tf, method = "average")
dend_tf <- as.dendrogram(hnorm_tf)

# 3) chromatin remodeler
dnorm_cr <- as.dist(1 - cor(expr_data_avg_nonImmune[rownames(expr_data_avg_nonImmune) %in% chromtainRM$V2, ], method = "pearson"))
hnorm_cr <- hclust(dnorm_cr, method = "average")
dend_cr <- as.dendrogram(hnorm_cr)

# 4) lncRNA
dnorm_lnc <- as.dist(1 - cor(expr_data_avg_nonImmune[rownames(expr_data_avg_nonImmune) %in% lncRNA$V1, ], method = "pearson"))
hnorm_lnc <- hclust(dnorm_lnc, method = "average")
dend_lnc <- as.dendrogram(hnorm_lnc)

# 5) cell surface
dnorm_cs <- as.dist(1 - cor(expr_data_avg_nonImmune[rownames(expr_data_avg_nonImmune) %in% cellSurface$V2, ], method = "pearson"))
hnorm_cs <- hclust(dnorm_cs, method = "average")
dend_cs <- as.dendrogram(hnorm_cs)

# N) module map
# module_map <- read.table(paste0("03-expression/merged/geneModule/", samplingPos, "/module_map.txt"), header = T, sep = "\t", row.names = 1, check.names = F)
# module_map_norm <- sweep(module_map, 1, rowSums(module_map, na.rm = T), "/")
# dnorm_mdm <- dist(module_map_norm)
# hnorm_mdm <- hclust(dnorm_mdm, method = "ward.D2", members = NULL)
# plot(as.dendrogram(hnorm_mdm), cex=.5, horiz=T)

pdf(paste0(OUT, "/hc_nonImmune.pdf"), width = 3, height = 14)

plot_hc(hnorm_all, dend_all, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(3, 90))
plot_hc(hnorm_tf, dend_tf, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(3, 90)); title(main = "TFs", line = 3)
plot_hc(hnorm_cr, dend_cr, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(3, 90)); title(main = "Chromatin remodelers", line = 3)
plot_hc(hnorm_lnc, dend_lnc, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(3, 90)); title(main = "LncRNAs", line = 3)
plot_hc(hnorm_cs, dend_cs, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(3, 90)); title(main = "Cell surface markers", line = 3)

dev.off()

# 1.3 Immune clgrp ----

# 1) all genes
dnorm_all <- as.dist(1 - cor(expr_data_avg_Immune, method = "pearson"))
hnorm_all <- hclust(dnorm_all, method = "average")
dend_all <- as.dendrogram(hnorm_all)

# 2) TF
dnorm_tf <- as.dist(1 - cor(expr_data_avg_Immune[rownames(expr_data_avg_Immune) %in% TF$V1, ], method = "pearson"))
hnorm_tf <- hclust(dnorm_tf, method = "average")
dend_tf <- as.dendrogram(hnorm_tf)

# 3) chromatin remodeler
dnorm_cr <- as.dist(1 - cor(expr_data_avg_Immune[rownames(expr_data_avg_Immune) %in% chromtainRM$V2, ], method = "pearson"))
hnorm_cr <- hclust(dnorm_cr, method = "average")
dend_cr <- as.dendrogram(hnorm_cr)

# 4) lncRNA
dnorm_lnc <- as.dist(1 - cor(expr_data_avg_Immune[rownames(expr_data_avg_Immune) %in% lncRNA$V1, ], method = "pearson"))
hnorm_lnc <- hclust(dnorm_lnc, method = "average")
dend_lnc <- as.dendrogram(hnorm_lnc)

# 5) cell surface
dnorm_cs <- as.dist(1 - cor(expr_data_avg_Immune[rownames(expr_data_avg_Immune) %in% chromtainRM$V2, ], method = "pearson"))
hnorm_cs <- hclust(dnorm_cs, method = "average")
dend_cs <- as.dendrogram(hnorm_cs)

# N) module map
# module_map <- read.table(paste0("03-expression/merged/geneModule/", samplingPos, "/module_map.txt"), header = T, sep = "\t", row.names = 1, check.names = F)
# module_map_norm <- sweep(module_map, 1, rowSums(module_map, na.rm = T), "/")
# dnorm_mdm <- dist(module_map_norm)
# hnorm_mdm <- hclust(dnorm_mdm, method = "ward.D2", members = NULL)
# plot(as.dendrogram(hnorm_mdm), cex=.5, horiz=T)

pdf(paste0(OUT, "/hc_Immune.pdf"), width = 3, height = 8)

plot_hc(hnorm_all, dend_all, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(1, 41))
plot_hc(hnorm_tf, dend_tf, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(1, 41))
plot_hc(hnorm_cr, dend_cr, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(1, 41))
plot_hc(hnorm_lnc, dend_lnc, cex.para = 0.8, mar.para = c(2,1.5,0,14) + 0.1, ylim.para = c(1, 41))
plot_hc(hnorm_cs, dend_cs, cex.para = 0.8, mar.para = c(2,1,1.5,14) + 0.1, ylim.para = c(1, 41))

dev.off()

# 3. TF expression ----
# bin expression
tf_bin <- apply(expr_data_tf_grp_norm, 2, function(x) {
  ### x <- expr_data_tf_grp_norm[, 20]
  bks <- seq(0, max(x), length.out = 9)
  y0 <- cut(x, breaks = bks, include.lowest = T, right = T)
  levels(y0) <- 1:length(levels(y0))
  yt <- table(y0)
  xxx <- data.frame(yt)
  xxx$y0 <- as.numeric(xxx$y0)
  xxx$group <- cumsum(c(T, xor(xxx$Freq[-1] == 0, xxx$Freq[-nrow(xxx)] == 0)))
  xxx <- subset(xxx, Freq != 0)
  if(length(unique(xxx$group)) == 1) {
    downhf <- floor(nrow(xxx) / 2) + 1
    xxx[downhf:nrow(xxx), "group"] <- xxx[downhf:nrow(xxx), "group"] + 0.5
  }
  
  xxx_LS <- lapply(split(xxx, xxx$group), function(x) {
    y <- x
    y$tg <- round(mean(y$y0))
    return(y)
  })
  xxx_DF <- do.call("rbind", xxx_LS)
  y1 <- xxx_DF[match(as.numeric(y0), xxx_DF$y0), "tg"]
  y2 <- bks[y1]
  out <- (y2 - min(y2)) / (max(y2) - min(y2))
  return(out)
})

tmp_t <- t(tf_bin)
tmp_t_sorted <- tmp_t[order(rowSums(tmp_t>0) < 10, apply(tmp_t, 1, which.max), -rowSums(tmp_t)), ]
tmp_n <- sum(rowSums(tmp_t>0) >= 10)

# Bin and scale according to gene expression distribution (from Nikos)
gene_cluster <- t(expr_data_tf_grp_norm)
gene_cluster_binned <- matrix(0,dim(gene_cluster)[1],dim(gene_cluster)[2])
for (i in 1:dim(gene_cluster)[1]) {
  x = rownames(gene_cluster)[i]
  y = as.matrix(gene_cluster[which(rownames(gene_cluster)%in%x), ])
  z = hist(y,n=10,plot=FALSE)
  #k[1,i]=table(z$density==0)["TRUE"]
  
  a=table(z$counts==0)["TRUE"]
  a[is.na(a)]=0
  b=length(z$breaks)
  
  if (a==0) {
    print("No gaps!")
    #gene_cluster_binned[i, ] <- 1
  }
  
  if (a==1) {
    l=which(z$counts==0)
    for(j in 1:dim(gene_cluster)[2]) {
      if (y[j] > z$breaks[l+1]) {
        gene_cluster_binned[i,j]=1
      }
    }
  }
  
  if (a>1) {
    for (j in 1:dim(gene_cluster)[2]) {
      for (k in 2:a) {
        n=which(z$counts==0)[k-1]
        l=which(z$counts==0)[k]
        if (y[j]<z$breaks[l+1] && y[j]>z$breaks[n+1]) {
          gene_cluster_binned[i,j]=(k-1)/a
        }
      }
      n=which(z$counts==0)[a]
      if (y[j]>z$breaks[n+1]) {
        gene_cluster_binned[i,j]=1
      }
    }
  }
  print(i)
}
rownames(gene_cluster_binned)=rownames(gene_cluster)
colnames(gene_cluster_binned)=colnames(gene_cluster)
#write.table(gene_cluster_binned,"cluster_gene_expression_digitized_n40.txt",quote=F,sep="\t")

gene_cluster_binned_sorted <- gene_cluster_binned[order(rowSums(gene_cluster_binned>0) < 10, apply(gene_cluster_binned, 1, which.max), -rowSums(gene_cluster_binned)), ]
###

# X. save ----
save.image(file = paste0(OUT, "/dendrogram.RData"))
