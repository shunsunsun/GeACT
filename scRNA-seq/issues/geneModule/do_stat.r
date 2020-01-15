
library("reshape2")
library("ggplot2")
library("cowplot")
library("ggalluvial")

# preprocess ----
fs <- list.files(path = "03-expression/merged/geneModule/", pattern = "geneModule.txt", recursive = T)
id <- gsub("/.*", "", fs)
fs_LS <- lapply(id, function(x) {
  y <- read.table(paste0("03-expression/merged/geneModule/", x , "/geneModule.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F)
  return(y)
})
names(fs_LS) <- id

module_table <- lapply(seq_along(fs_LS), function(i) {
  x <- fs_LS[[i]]
  y <- x[! duplicated(x$cluster), ]
  y$cond <- names(fs_LS)[i]
  y$cellNum <- as.numeric(sapply(strsplit(y$cond, split = "_"), "[", 1))
  y$seed <- as.numeric(sapply(strsplit(y$cond, split = "_"), "[", 2))
  y$method <- sapply(strsplit(y$cond, split = "_"), "[", 3)
  return(y)
})
module_table <- do.call("rbind", module_table)

# Module number / total genes / module size / avgCor ----
module_stat <- read.table("module_stat.txt", header = F, sep = "\t", stringsAsFactors = F)
dim(module_stat)
colnames(module_stat) <- c("cellNum", "seed", "method", "moduleNum", "geneNum")

pdf("module_stat.pdf", width = 5, height = 4)

ggplot(module_stat, aes(x = cellNum, y = moduleNum, color = factor(seed))) + geom_line(aes(linetype = method)) + geom_point() + 
  geom_vline(xintercept = 500, linetype = "dashed") + 
  scale_x_continuous(breaks = seq(100, 1700, by = 200)) + 
  scale_y_continuous(breaks = seq(0, 250, by = 50)) + 
  xlab("Cell number") + ylab("Module number") + labs(color = "Seed", linetype = "Clustering") + 
  theme(legend.position = c(0.5, 0.85), legend.spacing.y = unit(x = 0.01, units = "npc"), legend.box = "horizontal")

ggplot(module_stat, aes(x = cellNum, y = geneNum, color = factor(seed))) + geom_line(aes(linetype = method)) + geom_point() + 
  geom_vline(xintercept = 500, linetype = "dashed") + 
  scale_x_continuous(breaks = seq(100, 1700, by = 200)) + 
  scale_y_continuous(breaks = seq(0, 9000, by = 1500)) + 
  xlab("Cell number") + ylab("Total gene number in modules") + labs(color = "Seed", linetype = "Clustering") + 
  theme(legend.position = c(0.5, 0.65), legend.spacing.y = unit(x = 0.01, units = "npc"), legend.box = "horizontal")

ggplot(module_table, aes(x = factor(cellNum), y = size, fill = factor(seed))) + 
  geom_violin(show.legend = F) + facet_grid(method ~ ., scales = "free_y") + 
  xlab("Cell number") + ylab("Module size") + labs(fill = "Seed")

ggplot(module_table, aes(x = factor(cellNum), y = avgCor, fill = factor(seed))) + 
  geom_violin(show.legend = F) + facet_grid(method ~ ., scales = "fixed") + 
  xlab("Cell number") + ylab("Average correlation") + labs(fill = "Seed")

dev.off()

# Genes in modules overlapping ----

# single cond
do_geneOvlpScond <- function(fs_LS, num, seed1, seed2, hc_method = "average", show.legend = T) {
  g1 <- fs_LS[[paste(num, seed1, hc_method, sep = "_")]]$gene
  g2 <- fs_LS[[paste(num, seed2, hc_method, sep = "_")]]$gene
  n1 <- length(setdiff(g1, g2))
  n2 <- length(setdiff(g2, g1))
  bt <- length(intersect(g1, g2))
  ginmd_stat <-  data.frame(n1, n2, bt)
  colnames(ginmd_stat) <- c("OnlyIn1", "OnlyIn2", "Overlapped")
  ginmd_stat_melted <- melt(ginmd_stat, id.vars = NULL)
  ginmd_stat_melted$variable <- factor(ginmd_stat_melted$variable, levels = unique(ginmd_stat_melted$variable)[c(1,3,2)])
  print(ginmd_stat)
  print(colSums(ginmd_stat[, 1:3])/sum(ginmd_stat[, 1:3]))
  p <- ggplot(ginmd_stat_melted, aes(x = 1, y = value, fill = variable)) + geom_bar(stat = "identity", show.legend = show.legend) + 
    geom_text(aes(label = round(value/sum(value), 2)), position = position_stack(vjust = 0.5)) + 
    scale_y_continuous(expand = c(0, 0)) + ylab("Gene number") + labs(fill = NULL, title = paste(num, seed1, "~", seed2, hc_method)) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 4) + 
    coord_cartesian(clip = "off")
  return(p)
}

pdf("module_geneOvlpScond.pdf", width = 4, height = 4)
do_geneOvlpScond(fs_LS, num = 100, seed1 = 1, seed2 = 2, hc_method = "average")
do_geneOvlpScond(fs_LS, num = 500, seed1 = 1, seed2 = 2, hc_method = "average", show.legend = F)
do_geneOvlpScond(fs_LS, num = 1100, seed1 = 1, seed2 = 2, hc_method = "average", show.legend = F)
do_geneOvlpScond(fs_LS, num = 1700, seed1 = 1, seed2 = 2, hc_method = "average", show.legend = F)

do_geneOvlpScond(fs_LS, num = 100, seed1 = 1, seed2 = 2, hc_method = "complete")
do_geneOvlpScond(fs_LS, num = 500, seed1 = 1, seed2 = 2, hc_method = "complete", show.legend = F)
do_geneOvlpScond(fs_LS, num = 1100, seed1 = 1, seed2 = 2, hc_method = "complete", show.legend = F)
do_geneOvlpScond(fs_LS, num = 1700, seed1 = 1, seed2 = 2, hc_method = "complete", show.legend = F)
dev.off()

# multiple cond
do_geneOvlpMcond <- function(fs_LS, num1, num2, hc_method = "average") {
  ginmd_stat <- t(sapply(1:3, function(i) {
    g1 <- fs_LS[[paste(num1, i, hc_method, sep = "_")]]$gene
    g2 <- fs_LS[[paste(num2, i, hc_method, sep = "_")]]$gene
    n1 <- length(setdiff(g1, g2))
    n2 <- length(setdiff(g2, g1))
    bt <- length(intersect(g1, g2))
    return(c(n1, n2, bt))
  }))
  ginmd_stat <- as.data.frame(ginmd_stat)
  colnames(ginmd_stat) <- c("OnlyIn1", "OnlyIn2", "Overlapped")
  ginmd_stat$seed <- 1:3
  ginmd_stat_melted <- melt(ginmd_stat, id.vars = "seed")
  ginmd_stat_melted$variable <- factor(ginmd_stat_melted$variable, levels = unique(ginmd_stat_melted$variable)[c(1,3,2)])
  print(ginmd_stat)
  print(colSums(ginmd_stat[, 1:3])/sum(ginmd_stat[, 1:3]))
  ginmd_stat_ratio <- do.call("rbind", lapply(split(ginmd_stat_melted, ginmd_stat_melted$seed), function(x) { y <- x; y$ratio <- y$value / sum(y$value); return(y) }))
  p <- ggplot(ginmd_stat_melted, aes(x = seed, y = value, fill = variable)) + geom_bar(stat = "identity") + 
    geom_text(data = ginmd_stat_ratio, aes(label = round(ratio, 2)), position = position_stack(vjust = 0.5)) + 
    scale_y_continuous(expand = c(0, 0)) + xlab("Seed") + ylab("Gene number") + labs(fill = NULL, title = paste(num1, "~", num2, hc_method)) + 
    coord_cartesian(clip = "off")
  return(p)
}

pdf("module_geneOvlpMcond.pdf", width = 5, height = 5)
do_geneOvlpMcond(fs_LS, num1 = 100, num2 = 500, hc_method = "average")
do_geneOvlpMcond(fs_LS, num1 = 500, num2 = 1100, hc_method = "average")
do_geneOvlpMcond(fs_LS, num1 = 500, num2 = 1700, hc_method = "average")
do_geneOvlpMcond(fs_LS, num1 = 100, num2 = 500, hc_method = "complete")
do_geneOvlpMcond(fs_LS, num1 = 500, num2 = 1100, hc_method = "complete")
do_geneOvlpMcond(fs_LS, num1 = 500, num2 = 1700, hc_method = "complete")
dev.off()

# gene alluvial ----

# 500 vs 1100
t1 <- fs_LS[["500_1_average"]]
t1$cluster <- factor(t1$cluster, levels = unique(t1$cluster))
t2 <- fs_LS[["1100_1_average"]]
t2$cluster <- gsub("M", "S", t2$cluster)
t2$cluster <- factor(t2$cluster, levels = unique(t2$cluster))
tx <- merge(t1[, c(1:3,5,7)], t2[, c(1:3,5,7)], by = "gene", sort = F, all = T)
tx_stat <- as.data.frame(table(tx[, c("cluster.x", "cluster.y")], useNA = "ifany"))
tx_stat <- subset(tx_stat, Freq >= 1)
tx_stat <- tx_stat[order(tx_stat$Freq, decreasing = T), ]
rownames(tx_stat) <- NULL

# by 1
t1_rmdup <- t1[! duplicated(t1$cluster), ]
t1_rmdup <- t1_rmdup[order(t1_rmdup$avgCor, decreasing = T), ]
tx_stat_sub1 <- subset(tx_stat, cluster.x %in% t1_rmdup$cluster[1:10])
# by 2
t2_rmdup <- t2[! duplicated(t2$cluster), ]
t2_rmdup <- t2_rmdup[order(t2_rmdup$avgCor, decreasing = T), ]
tx_stat_sub2 <- subset(tx_stat, cluster.y %in% t2_rmdup$cluster[1:10])

pdf("module_alluvial_500_1100.pdf", width = 5, height = 10)
ggplot(tx_stat, aes(axis1 = cluster.x, axis2 = cluster.y, y = Freq)) + geom_alluvium(aes(fill = cluster.x), width = 1/12, show.legend = F) + 
  geom_stratum(fill = "skyblue", color = "grey90", width = 1/12) + 
  scale_x_discrete(limits = c("cluster.x", "cluster.y"), labels = c("N500", "N1100"), expand = c(0.05, 0.05)) + ylab("Gene number") + 
  theme(aspect.ratio = 2)

ggplot(tx_stat_sub1, aes(axis1 = cluster.x, axis2 = cluster.y, y = Freq)) + geom_alluvium(aes(fill = cluster.x), width = 1/12, show.legend = F) + 
  geom_stratum(fill = "skyblue", color = "grey90", width = 1/12) + 
  geom_label(stat = "stratum", infer.label = T) + 
  scale_x_discrete(limits = c("cluster.x", "cluster.y"), labels = c("N500", "N1100"), expand = c(0.05, 0.05)) + ylab("Gene number") + 
  theme(aspect.ratio = 2)

ggplot(tx_stat_sub2, aes(axis1 = cluster.x, axis2 = cluster.y, y = Freq)) + geom_alluvium(aes(fill = cluster.y), width = 1/12, show.legend = F) + 
  geom_stratum(fill = "skyblue", color = "grey90", width = 1/12) + 
  geom_label(stat = "stratum", infer.label = T) + 
  scale_x_discrete(limits = c("cluster.x", "cluster.y"), labels = c("N500", "N1100"), expand = c(0.05, 0.05)) + ylab("Gene number") + 
  theme(aspect.ratio = 2)
dev.off()

# 500 vs 1700
t1 <- fs_LS[["500_1_average"]]
t1$cluster <- factor(t1$cluster, levels = unique(t1$cluster))
t2 <- fs_LS[["1700_1_average"]]
t2$cluster <- gsub("M", "S", t2$cluster)
t2$cluster <- factor(t2$cluster, levels = unique(t2$cluster))
tx <- merge(t1[, c(1:3,5,7)], t2[, c(1:3,5,7)], by = "gene", sort = F, all = T)
tx_stat <- as.data.frame(table(tx[, c("cluster.x", "cluster.y")], useNA = "ifany"))
tx_stat <- subset(tx_stat, Freq >= 1)
tx_stat <- tx_stat[order(tx_stat$Freq, decreasing = T), ]
rownames(tx_stat) <- NULL

# by 1
t1_rmdup <- t1[! duplicated(t1$cluster), ]
t1_rmdup <- t1_rmdup[order(t1_rmdup$avgCor, decreasing = T), ]
tx_stat_sub1 <- subset(tx_stat, cluster.x %in% t1_rmdup$cluster[1:10])
# by 2
t2_rmdup <- t2[! duplicated(t2$cluster), ]
t2_rmdup <- t2_rmdup[order(t2_rmdup$avgCor, decreasing = T), ]
tx_stat_sub2 <- subset(tx_stat, cluster.y %in% t2_rmdup$cluster[1:10])

pdf("module_alluvial_500_1700.pdf", width = 5, height = 10)
ggplot(tx_stat, aes(axis1 = cluster.x, axis2 = cluster.y, y = Freq)) + geom_alluvium(aes(fill = cluster.x), width = 1/12, show.legend = F) + 
  geom_stratum(fill = "skyblue", color = "grey90", width = 1/12) + 
  scale_x_discrete(limits = c("cluster.x", "cluster.y"), labels = c("N500", "N1700"), expand = c(0.05, 0.05)) + ylab("Gene number") + 
  theme(aspect.ratio = 2)

ggplot(tx_stat_sub1, aes(axis1 = cluster.x, axis2 = cluster.y, y = Freq)) + geom_alluvium(aes(fill = cluster.x), width = 1/12, show.legend = F) + 
  geom_stratum(fill = "skyblue", color = "grey90", width = 1/12) + 
  geom_label(stat = "stratum", infer.label = T) + 
  scale_x_discrete(limits = c("cluster.x", "cluster.y"), labels = c("N500", "N1700"), expand = c(0.05, 0.05)) + ylab("Gene number") + 
  theme(aspect.ratio = 2)

ggplot(tx_stat_sub2, aes(axis1 = cluster.x, axis2 = cluster.y, y = Freq)) + geom_alluvium(aes(fill = cluster.y), width = 1/12, show.legend = F) + 
  geom_stratum(fill = "skyblue", color = "grey90", width = 1/12) + 
  geom_label(stat = "stratum", infer.label = T) + 
  scale_x_discrete(limits = c("cluster.x", "cluster.y"), labels = c("N500", "N1700"), expand = c(0.05, 0.05)) + ylab("Gene number") + 
  theme(aspect.ratio = 2)
dev.off()
