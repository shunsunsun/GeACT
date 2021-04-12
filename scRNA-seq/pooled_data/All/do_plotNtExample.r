# plot network example
setwd("~/lustre/06-Human_cell_atlas/pooled_data/All/")

library("ggplot2")
library("cowplot")
suppressMessages(library("igraph"))

edges <- data.frame(TF = "FOXL1", target = c("ETS2", "F2R", "PDLIM1", "EHD3"), stringsAsFactors = F)
nodes <- data.frame(id = unique(c(edges$TF, edges$target)), stringsAsFactors = F)

nodes$color <- c("limegreen", "tomato", "orange", "royalblue1", "mediumpurple1")
nodes$label <- NA
nodes$frame.color <- "grey90"

g <- graph_from_data_frame(edges, directed = T, vertices = nodes)
list.vertex.attributes(g)
list.edge.attributes(g)
coords <- layout_as_tree(g, root = "FOXL1", "all")
rownames(coords) <- as_ids(V(g))

coords[1, ] <- c(2.5, 2)
coords[-1, 1] <- 1:4
coords[-1, 2] <- 1

pdf("03-expression/merged/geneModule/network_simple.pdf", width = 4, height = 4)
par(mar=c(0,0,0,0))
plot.igraph(g, layout = coords, asp = 1, 
            vertex.size = 35, 
            vertex.shape = "circle", 
            #vertex.frame.color = "white", 
            vertex.label.color = "black", 
            #vertex.label.cex = 0.2, 
            edge.color = "black", 
            edge.curved = 0, 
            #edge.arrow.size= 0.2, 
            xlim = c(-1, 1), ylim = c(-1, 1)
)
dev.off()
