
do_plotNt <- function(gs_case, col.low = "dark blue", col.mid = "white", col.high = "dark red", show.all.nodes = F, print.edges = F, print.nodes = F, do.plot = T) {
  edges <- subset(network_ftd_DF, geneSet %in% paste(gs_case, c("up", "down"), sep = "."), c("TF", "target"))
  nodes <- data.frame(gene = unique(c(edges$TF, edges$target)), stringsAsFactors = F)
  rownames(nodes) <- nodes$gene
  expr.markers_ftd_sub <- subset(expr.markers_ftd, cluster == gs_case)
  nodes$avg_logFC <- expr.markers_ftd_sub$avg_logFC[match(nodes$gene, expr.markers_ftd_sub$gene)]
  nodes <- nodes[order(nodes$avg_logFC), ]
  #range(nodes$avg_logFC, na.rm = T)
  tf_id <- unique(edges$TF)
  tg_id <- setdiff(nodes$gene, tf_id)
  if(length(tg_id) == 0) {
    cat("No target genes.\n")
    return(NULL)
  }
  nodes$color <- colorRampPalette(c(col.low, col.mid, col.high))(100)[cut(nodes$avg_logFC, breaks = seq(-2.5, 2.5, length.out = 101))]
  nodes$color[nodes$avg_logFC < -2.5] <- col.low
  nodes$color[nodes$avg_logFC > 2.5] <- col.high
  nodes$color[is.na(nodes$avg_logFC)] <- "grey80"
  nodes$label <- nodes$gene
  if(! show.all.nodes) {
    nodes$label[(! nodes$label %in% tf_id) & nodes$avg_logFC < log(3) & nodes$avg_logFC > -log(3)] <- NA
  }
  nodes$label.dist <- ifelse(nodes$label %in% tf_id, 0, 2.5)
  nodes$frame.color <- "grey90"
  
  if(print.edges) { print(edges) }
  if(print.nodes) { print(nodes) }
  
  g <- graph_from_data_frame(edges, directed = T, vertices = nodes)
  #list.vertex.attributes(g)
  #list.edge.attributes(g)
  
  set.seed(1234)
  coords <- layout_with_fr(g)
  rownames(coords) <- as_ids(V(g))
  coords[tf_id, 1] <- cos(seq(0,2*pi*(1 - 1 / length(tf_id)), length.out = length(tf_id))) * 0.2
  coords[tf_id, 2] <- sin(seq(0,2*pi*(1 - 1 / length(tf_id)), length.out = length(tf_id))) * 0.2
  coords[tg_id, 1] <- cos(seq(0,2*pi*(1 - 1 / length(tg_id)), length.out = length(tg_id))) * 1
  coords[tg_id, 2] <- sin(seq(0,2*pi*(1 - 1 / length(tg_id)), length.out = length(tg_id))) * 1
  
  #return(list(g = g, coords = coords))
  
  if(do.plot) {
    if(is.null(g)) {
      cat("No network to be shown.\n")
    } else {
      par(mar=c(0,0,2,0))
      plot.igraph(g, layout = coords, asp = 1, 
                  vertex.size = sqrt(degree(g)) * 2.5, #vertex.size = 6, 
                  vertex.shape = "circle", 
                  #vertex.frame.color = "white", 
                  vertex.label.color = "black", 
                  #vertex.label.cex = 0.2, 
                  vertex.label.degree = 0, 
                  edge.curved = 0, 
                  edge.arrow.size= 0.4, 
                  xlim = c(-1, 1.2), ylim = c(-1.2, 0.95), 
                  main = gs_case
      )
      fields::image.plot(col = colorRampPalette(c(col.low, col.mid, col.high))(100), legend.only = T, 
                         legend.line = -2.05, legend.cex = 0.8, legend.width = 0.8, legend.lab = "Log (Fold change)", 
                         axis.args = list(mgp= c(3,0,0), cex.axis = 0.7, tck = -0.2),
                         zlim = c(-2.5, 2.5), horizontal = T, legend.shrink = 0.1, smallplot=c(.15, .25, 0.05, 0.075))
      par(mar=c(5,4,4,2) + 0.1)
    }
  }
}


do_plotNt_alt <- function(gs_case, type = "both", col.low = "dark blue", col.mid = "white", col.high = "dark red", tf.pos.r = 0.2, 
                          xlim = c(-1, 1), ylim = c(-1, 1), edge.color = "grey80", color_TF = T, 
                          show.node.names = "TF_DEG", print.edges = F, print.nodes = F, do.plot = T) {
  if(type == "both") {
    gs_case.used <- paste(gs_case, c("up", "down"), sep = ".")
  } else if (type == "up") {
    gs_case.used <- paste(gs_case, c("up"), sep = ".")
  } else if (type == "down") {
    gs_case.used <- paste(gs_case, c("down"), sep = ".")
  } else {
    stop("type should be: both|up|down")
  }
  edges <- subset(network_ftd_DF, geneSet %in% gs_case.used, c("TF", "target", "geneSet"))
  nodes <- data.frame(gene = unique(c(edges$TF, edges$target)), stringsAsFactors = F)
  rownames(nodes) <- nodes$gene
  expr.markers_ftd_sub <- subset(expr.markers_ftd, cluster == gs_case)
  nodes$avg_logFC <- expr.markers_ftd_sub$avg_logFC[match(nodes$gene, expr.markers_ftd_sub$gene)]
  nodes <- nodes[order(nodes$avg_logFC), ]
  #range(nodes$avg_logFC, na.rm = T)
  tf_id <- unique(edges$TF)
  tg_id <- setdiff(nodes$gene, tf_id)
  if(length(tg_id) == 0) {
    cat("No target genes.\n")
    return(NULL)
  }
  nodes$color <- colorRampPalette(c(col.low, col.mid, col.high))(100)[cut(nodes$avg_logFC, breaks = seq(-2.5, 2.5, length.out = 101))]
  nodes$color[nodes$avg_logFC < -2.5] <- col.low
  nodes$color[nodes$avg_logFC > 2.5] <- col.high
  nodes$color[is.na(nodes$avg_logFC)] <- "grey80"
  nodes$frame.color <- "grey90"
  if(color_TF) {
    TF_up <- unique(subset(edges, grepl("up$", geneSet), "TF", drop = T))
    TF_down <- unique(subset(edges, grepl("down$", geneSet), "TF", drop = T))
    TF_both <- c(TF_up, TF_down)[duplicated(c(TF_up, TF_down))]
    nodes$frame.color[nodes$gene %in% TF_up] <- "red"
    nodes$frame.color[nodes$gene %in% TF_down] <- "blue"
    nodes$frame.color[nodes$gene %in% TF_both] <- "purple"
  }
  nodes$label <- nodes$gene
  if(show.node.names == "TF") {
    nodes$label[! nodes$label %in% tf_id] <- NA
  } else if(show.node.names == "TF_DEG") {
    nodes$label[(! nodes$label %in% tf_id) & nodes$avg_logFC < log(3) & nodes$avg_logFC > -log(3)] <- NA
  }
  nodes$label.dist <- ifelse(nodes$label %in% tf_id, 0, 2.5)
  
  if(print.edges) { print(edges) }
  if(print.nodes) { print(nodes) }
  
  g <- graph_from_data_frame(edges, directed = T, vertices = nodes)
  #list.vertex.attributes(g)
  #list.edge.attributes(g)
  
  set.seed(1234)
  coords <- layout_with_fr(g)
  rownames(coords) <- as_ids(V(g))
  coords[tf_id, 1] <- cos(seq(0,2*pi*(1 - 1 / length(tf_id)), length.out = length(tf_id))) * tf.pos.r
  coords[tf_id, 2] <- sin(seq(0,2*pi*(1 - 1 / length(tf_id)), length.out = length(tf_id))) * tf.pos.r
  coords[tg_id, 1] <- cos(seq(0,2*pi*(1 - 1 / length(tg_id)), length.out = length(tg_id))) * 1
  coords[tg_id, 2] <- sin(seq(0,2*pi*(1 - 1 / length(tg_id)), length.out = length(tg_id))) * 1
  
  #return(list(g = g, coords = coords))
  
  if(do.plot) {
    if(is.null(g)) {
      cat("No network to be shown.\n")
    } else {
      par(mar=c(0,0,2,0))
      plot.igraph(g, layout = coords, asp = 1, 
                  vertex.size = sqrt(degree(g)) * 2.5, #vertex.size = 6, 
                  vertex.shape = "circle", 
                  #vertex.frame.color = "white", 
                  vertex.label.color = "black", 
                  #vertex.label.cex = 0.2, 
                  vertex.label.degree = 0, 
                  edge.color = edge.color, 
                  edge.curved = 0, 
                  edge.arrow.size= 0.4, 
                  xlim = xlim, ylim = ylim
      )
      fields::image.plot(col = colorRampPalette(c(col.low, col.mid, col.high))(100), legend.only = T, 
                         legend.line = -2.1, legend.cex = 1, legend.width = 0.8, legend.lab = "Log (Fold change)", 
                         axis.args = list(at = c(-2, 0, 2), mgp= c(3, 0.2, 0), cex.axis = 1, tck = -0.2),
                         zlim = c(-2.5, 2.5), horizontal = T, legend.shrink = 0.1, smallplot=c(.4, .6, 0.05, 0.075))
      par(mar=c(5,4,4,2) + 0.1)
    }
  }
}

do_shellPlot <- function(gs_case, type = "both", col.low = "dark blue", col.mid = "white", col.high = "dark red", tf.pos.r = 0.2, y.dist = 0.1, 
                         xlim = c(-1, 1), ylim = c(-1, 1), edge.color = "grey80", color_TF = T, 
                         show.node.names = "TF_DEG", print.edges = F, print.nodes = F, do.plot = T) {
  if(type == "both") {
    gs_case.used <- paste(gs_case, c("up", "down"), sep = ".")
  } else if (type == "up") {
    gs_case.used <- paste(gs_case, c("up"), sep = ".")
  } else if (type == "down") {
    gs_case.used <- paste(gs_case, c("down"), sep = ".")
  } else {
    stop("type should be: both|up|down")
  }
  edges <- subset(network_ftd_DF, geneSet %in% gs_case.used, c("TF", "target", "geneSet"))
  nodes <- data.frame(gene = unique(c(edges$TF, edges$target)), stringsAsFactors = F)
  rownames(nodes) <- nodes$gene
  expr.markers_ftd_sub <- subset(expr.markers_ftd, cluster == gs_case)
  nodes$avg_logFC <- expr.markers_ftd_sub$avg_logFC[match(nodes$gene, expr.markers_ftd_sub$gene)]
  nodes <- nodes[order(nodes$avg_logFC), ]
  #range(nodes$avg_logFC, na.rm = T)
  tf_id <- unique(edges$TF)
  tg_id <- setdiff(nodes$gene, tf_id)
  if(length(tg_id) == 0) {
    cat("No target genes.\n")
    return(NULL)
  }
  nodes$color <- colorRampPalette(c(col.low, col.mid, col.high))(100)[cut(nodes$avg_logFC, breaks = seq(-2.5, 2.5, length.out = 101))]
  nodes$color[nodes$avg_logFC < -2.5] <- col.low
  nodes$color[nodes$avg_logFC > 2.5] <- col.high
  nodes$color[is.na(nodes$avg_logFC)] <- "grey80"
  nodes$frame.color <- "grey90"
  if(color_TF) {
    TF_up <- unique(subset(edges, grepl("up$", geneSet), "TF", drop = T))
    TF_down <- unique(subset(edges, grepl("down$", geneSet), "TF", drop = T))
    TF_both <- c(TF_up, TF_down)[duplicated(c(TF_up, TF_down))]
    nodes$frame.color[nodes$gene %in% TF_up] <- "red"
    nodes$frame.color[nodes$gene %in% TF_down] <- "blue"
    nodes$frame.color[nodes$gene %in% TF_both] <- "purple"
  }
  nodes$label <- nodes$gene
  if(show.node.names == "TF") {
    nodes$label[! nodes$label %in% tf_id] <- NA
  } else if(show.node.names == "TF_DEG") {
    nodes$label[(! nodes$label %in% tf_id) & nodes$avg_logFC < log(3) & nodes$avg_logFC > -log(3)] <- NA
  }
  nodes$label.dist <- ifelse(nodes$label %in% tf_id, 0, 2.5)
  
  if(print.edges) { print(edges) }
  if(print.nodes) { print(nodes) }
  
  g <- graph_from_data_frame(edges, directed = T, vertices = nodes)
  #list.vertex.attributes(g)
  #list.edge.attributes(g)
  
  set.seed(1234)
  coords <- layout_with_fr(g)
  rownames(coords) <- as_ids(V(g))
  # split
  tg_id_up <- intersect(tg_id, nodes$gene[nodes$avg_logFC > 0])
  tf_id_up <- unique(subset(edges, target %in% tg_id_up, "TF", drop = T))
  tg_id_dw <- intersect(tg_id, nodes$gene[nodes$avg_logFC < 0])
  tf_id_dw <- unique(subset(edges, target %in% tg_id_dw, "TF", drop = T))
  tf_id_bt <- c(tf_id_up, tf_id_dw)[duplicated(c(tf_id_up, tf_id_dw))]
  tf_id_up <- setdiff(tf_id_up, tf_id_bt)
  tf_id_dw <- setdiff(tf_id_dw, tf_id_bt)
  # pos (up)
  coords[tf_id_up, 1] <- cos(seq(0, 1*pi, length.out = length(tf_id_up) + 2)[- c(1, length(tf_id_up) + 2)]) * tf.pos.r * 1.25
  coords[tf_id_up, 2] <- sin(seq(0, 1*pi, length.out = length(tf_id_up) + 2)[- c(1, length(tf_id_up) + 2)]) * tf.pos.r + y.dist / 2
  coords[tg_id_up, 1] <- cos(seq(0, 1*pi, length.out = length(tg_id_up) + 2)[- c(1, length(tg_id_up) + 2)]) * 1
  coords[tg_id_up, 2] <- sin(seq(0, 1*pi, length.out = length(tg_id_up) + 2)[- c(1, length(tg_id_up) + 2)]) * 1 + y.dist / 2
  # pos (dw)
  coords[tf_id_dw, 1] <- cos(-seq(0, 1*pi, length.out = length(tf_id_dw) + 2)[- c(1, length(tf_id_dw) + 2)]) * tf.pos.r * 1.25
  coords[tf_id_dw, 2] <- sin(-seq(0, 1*pi, length.out = length(tf_id_dw) + 2)[- c(1, length(tf_id_dw) + 2)]) * tf.pos.r - y.dist / 2
  coords[tg_id_dw, 1] <- cos(-seq(0, 1*pi, length.out = length(tg_id_dw) + 2)[- c(1, length(tg_id_dw) + 2)]) * 1
  coords[tg_id_dw, 2] <- sin(-seq(0, 1*pi, length.out = length(tg_id_dw) + 2)[- c(1, length(tg_id_dw) + 2)]) * 1 - y.dist / 2
  # pos (both)
  if(length(tf_id_bt) > 0) {
    coords[tf_id_bt, 1] <- cos(seq(0,2*pi*(1 - 1 / length(tf_id_bt)), length.out = length(tf_id_bt))) * 0.1
    coords[tf_id_bt, 2] <- sin(seq(0,2*pi*(1 - 1 / length(tf_id_bt)), length.out = length(tf_id_bt))) * 0.1
  }
  
  #return(list(g = g, coords = coords))
  
  if(do.plot) {
    if(is.null(g)) {
      cat("No network to be shown.\n")
    } else {
      par(mar=c(0,0,0,0))
      plot.igraph(g, layout = coords, asp = 1, 
                  vertex.size = sqrt(degree(g)) * 2.5, #vertex.size = 6, 
                  vertex.shape = "circle", 
                  #vertex.frame.color = "white", 
                  vertex.label.color = "black", 
                  #vertex.label.cex = 0.2, 
                  vertex.label.degree = 0, 
                  edge.color = edge.color, 
                  edge.curved = 0, 
                  edge.arrow.size= 0.4, 
                  xlim = xlim, ylim = ylim
      )
      fields::image.plot(col = colorRampPalette(c(col.low, col.mid, col.high))(100), legend.only = T, 
                         legend.line = -2.1, legend.cex = 1, legend.width = 0.8, legend.lab = "Log (Fold change)", 
                         axis.args = list(at = c(-2, 0, 2), mgp= c(3, 0.2, 0), cex.axis = 1, tck = -0.2),
                         zlim = c(-2.5, 2.5), horizontal = T, legend.shrink = 0.1, smallplot=c(.4, .6, 0.05, 0.075))
      par(mar=c(5,4,4,2) + 0.1)
    }
  }
}
