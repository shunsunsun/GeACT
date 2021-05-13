#!/usr/bin/Rscript

if(length(commandArgs(trailingOnly = T)) != 3) {
  cat("Usage:   Rscript do_exprPlot.r tissue dimension_type gene\n")
  cat("Example: Rscript do_exprPlot.r small_intestine tSNE COL1A2\n")
  q("no")
}

tissue <- commandArgs(trailingOnly = T)[1]
dim_type <- commandArgs(trailingOnly = T)[2]
gene <- commandArgs(trailingOnly = T)[3]

# check gene ID
gene_ID <- read.table(file = "data/gene_ID_mapping.txt", header = T, sep = "\t", stringsAsFactors = F)
if(grepl("^ENSG", gene)) {
  if(grepl("\\.[0-9]+$", gene)) {
    gene_all <- gene_ID$ensembl
  } else {
    gene_all <- gene_ID$ensembl_short
  }
} else {
  gene_all <- gene_ID$symbol
}
if(! gene %in% gene_all) {
  stop("Unrecognized gene ID.")
}
gene <- gene_ID$symbol[gene_all == gene]

OUT <- "output"

outfile <- paste0(OUT, "/", tissue, "/", dim_type, "/", gene, ".png")
if(! file.exists(outfile)) {
  suppressMessages(library("arrow"))
  library("ggplot2")
  suppressMessages(library("cowplot"))
  suppressMessages(library("R.devices"))
  #suppressMessages(library("plotly"))
  
  cellMatr <- read_feather(file = paste0("data/", tissue, "_cellMatr.feather"))
  cellMatr <- as.data.frame(cellMatr)
  rownames(cellMatr) <- read.table(file = paste0("data/", tissue, "_cellMatr.cell"), header = F, sep = "\t", stringsAsFactors = F)$V1
  
  gp <- ggplot(cellMatr, aes_string(x = paste0(dim_type, "_1_ali"), y = paste0(dim_type, "_2_ali"))) + geom_point(aes_string(color = gene), alpha = 0.9, size = 0.9) + 
    scale_color_gradient(low = "grey", high = "limegreen", limits = c(0, max(c(cellMatr[, gene], 5)))) + 
    xlab(NULL) + ylab(NULL) + labs(color = parse(text = "Log[2]~(CPM+1)")) + ggtitle(gene) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + 
    theme(aspect.ratio = 1) + 
    theme(plot.title = element_text(face = "italic"))
  
  dir.create(path = paste0(OUT, "/", tissue, "/", dim_type), showWarnings = F, recursive = T)
  invisible(suppressGraphics(ggsave(filename = outfile, plot = gp, device = "png", width = 8, height = 6, units = "in")))
}
#cat("Output file:", outfile, "\n")

#gl <- ggplotly(gp)
#htmlwidgets::saveWidget(widget = gl, file = paste0(tissue, "_", gene, "_", dim_type, ".html"))
