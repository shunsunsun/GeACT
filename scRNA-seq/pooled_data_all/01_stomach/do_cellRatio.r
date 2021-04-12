# cell ratio analysis
setwd("~/lustre/06-Human_cell_atlas/pooled_data_all/01_stomach/")

library("reshape2")
library("RColorBrewer")
library("ggplot2")
suppressMessages(library("cowplot"))

samplingPos <- "."
OUT <- paste0("03-expression/merged/cellRatio/", samplingPos)
dir.create(OUT, showWarnings = F, recursive = T)

#load(file = paste0(OUT, "/cellRatio.RData"))

# 1. preprocess ----
# cell meta
cellMetaData <- read.table(file = "cell_metatable_filtered_plus.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
cellMetaData[cellMetaData$group %in% c("B", "DC/Macrophage", "T"), "group"] <- "Immune"
cellMetaData$samplingPos <- Hmisc::capitalize(cellMetaData$samplingPos)

# cell type meta
ctMetaData <- read.table(file = "cellType_metatable.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
ctMetaData <- subset(ctMetaData, stage == "14w")[-4]
ctMetaData <- merge(ctMetaData, unique(cellMetaData[, c("ident", "group")]), by = "ident", sort = F)

cellMetaData$samplingPos <- factor(cellMetaData$samplingPos, levels = c("Fundus", "Body", "Antrum"))
cellMetaData$ident <- factor(cellMetaData$ident, levels = ctMetaData$ident)
cellMetaData$group <- factor(cellMetaData$group, levels = unique(ctMetaData$group))

# color
load("../../pooled_data/All/03-expression/merged/cellCluster/ct_color.RData")
cg_color <- c("tomato", "orange", "#FB9A99", "#00BF74", "#8258FA",
              "#D0A9F5", "purple2", "#D358F7", "#419FDE", "#F2E100",
              ct_color["Ovary"], ct_color["Testis"], "firebrick3", "grey80", 
              "#8258FA", "skyblue", "#B3B3B3")
names(cg_color) <- c("Epithelial", "Endothelial", "Smooth muscle", "Fibroblast", "B", 
                     "DC/Macrophage", "NKT", "T", "Glial", "FGC", "Granulosa", "Sertoli", "Erythrocyte", "Other", 
                     "Immune", "CACNA1A", "Unknown")
cg_color <- cg_color[levels(cellMetaData$group)]

# 2. all cells for analysis ----
cellRatio_LS <- lapply(unique(cellMetaData$samplingPos), function(i) {
  x <- subset(cellMetaData, samplingPos == i)
  y0 <- table(x[, c("stage", "group")])
  #y1 <- sweep(y0, MARGIN = 1, STATS = rowSums(y0), "/")
  y1 <- y0
  y <- melt(y1)
  y$samplingPos <- i
  return(y)
})
cellRatio_DF <- do.call("rbind", cellRatio_LS)

pdf(paste0(OUT, "/cellRatio_all.pdf"), width = 6, height = 4)

gp <- ggplot(cellRatio_DF, aes(x = stage, y = value, fill = group)) + geom_bar(stat = "identity", position = "fill") + 
  facet_grid(. ~ samplingPos) + 
  scale_fill_manual(values = cg_color) + scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) + 
  labs(fill = NULL) + xlab(NULL) + ylab("Cell ratio to all cells") + 
  theme(strip.background = element_blank())
pg <- ggplotGrob(gp)
for(i in which(grepl("strip-t", pg$layout$name))){
  pg$grobs[[i]]$layout$clip <- "off"
}
grid::grid.draw(pg)

dev.off()

# 3. fibroblast ----
cellMetaData_sub <- subset(cellMetaData, group == "Fibroblast")
cellMetaData_sub$ident <- factor(cellMetaData_sub$ident)
cellRatio_sub_LS <- lapply(unique(cellMetaData_sub$samplingPos), function(i) {
  x <- subset(cellMetaData_sub, samplingPos == i)
  y0 <- table(x[, c("stage", "ident")])
  #y1 <- sweep(y0, MARGIN = 1, STATS = rowSums(y0), "/")
  y1 <- y0
  y <- melt(y1)
  y$samplingPos <- i
  return(y)
})
cellRatio_sub_DF <- do.call("rbind", cellRatio_sub_LS)

pdf(paste0(OUT, "/cellRatio_fibro.pdf"), width = 6, height = 4)

gp <- ggplot(cellRatio_sub_DF, aes(x = stage, y = value, fill = ident)) + geom_bar(stat = "identity", position = "fill") + 
  facet_grid(. ~ samplingPos) + 
  scale_fill_manual(values = c(colorRampPalette(brewer.pal(8, "Set2"))(8)[-8], "slateblue1", "palegreen2")) + scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) + 
  labs(fill = NULL) + xlab(NULL) + ylab("Cell ratio to all fibroblasts") + 
  theme(strip.background = element_blank())
pg <- ggplotGrob(gp)
for(i in which(grepl("strip-t", pg$layout$name))){
  pg$grobs[[i]]$layout$clip <- "off"
}
grid::grid.draw(pg)

dev.off()

save.image(file = paste0(OUT, "/cellRatio.RData"))
