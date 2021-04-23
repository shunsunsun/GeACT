
do_plotTrack <- function(case_gene, highlight.dsl = T, highlight.dsl.width = 0.01, vmin = -0.6, vmax = 0.6, margin.right = 30) {
  #cat(">", case_gene, "\n")
  dsl_case <- subset(dsl_s, gene == case_gene)
  if(nrow(dsl_case) == 0) {
    cat("No differential significant links.\n")
    next
  }
  dsl_case <- dsl_case[order(dsl_case$peakLeft), ]
  tss_pos <- subset(gene_pos, gene == case_gene, "tss", drop = T)
  
  link_case_11_14 <- subset(link_11_14, gene == case_gene)
  link_case_11_14$sig <- link_case_11_14$peak %in% subset(link_11_14_ftd, gene == case_gene, "peak", drop = T)
  link_case_19_22 <- subset(link_19_22, gene == case_gene)
  link_case_19_22$sig <- link_case_19_22$peak %in% subset(link_19_22_ftd, gene == case_gene, "peak", drop = T)
  link_case <- rbind(link_case_11_14, link_case_19_22)
  link_case$stage <- rep(c("11-14w", "19-22w"), c(nrow(link_case_11_14), nrow(link_case_19_22)))
  link_case$peakChr <- gsub(":.*", "", link_case$peak)
  link_case$peakLeft <- as.numeric(gsub("-.*", "", gsub(".*:", "", link_case$peak)))
  link_case$peakRight <- as.numeric(gsub(".*-", "", gsub(".*:", "", link_case$peak)))
  link_case$peakCenter <- (link_case$peakLeft + link_case$peakRight) / 2
  link_case$tss_pos <- tss_pos
  link_case$direction <- ifelse(link_case$peakCenter > link_case$tss_pos, 1, -1)
  #
  link_case$Correlation <- pmin(link_case$Correlation, vmax)
  link_case$Correlation <- pmax(link_case$Correlation, vmin)
  #
  xlim.min <- min(c(link_case$peakLeft, link_case$peakRight, tss_pos))
  xlim.max <- max(c(link_case$peakLeft, link_case$peakRight, tss_pos))
  # make gp
  gp0 <- ggplot(link_case, aes(x = peakCenter, xend = tss_pos, y = 0, yend = 0, color = Correlation)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    coord_cartesian(xlim = c(xlim.min, xlim.max), ylim = c(0, 1), clip = "off") + 
    xlab(NULL) + ylab(NULL) + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = margin(7, margin.right, 7, 7, unit = "pt")) + 
    theme(plot.title = element_text(face = "plain")) + 
    scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = c(vmin, vmax)) + 
    theme(legend.position = "none")
  # 11-14w
  gp1 <- gp0 + ylab("11-14w") + ggtitle(case_gene)
  link_case_sub <- subset(link_case, stage == "11-14w" & (! sig) & direction == -1)
  if(nrow(link_case_sub) > 0) {
    gp1 <- gp1 + geom_curve(data = link_case_sub, curvature = -0.5, ncp = 10, color = "grey90")
  }
  link_case_sub <- subset(link_case, stage == "11-14w" & (! sig) & direction == 1)
  if(nrow(link_case_sub) > 0) {
    gp1 <- gp1 + geom_curve(data = link_case_sub, curvature = 0.5, ncp = 10, color = "grey90")
  }
  link_case_sub <- subset(link_case, stage == "11-14w" & sig & direction == -1)
  if(nrow(link_case_sub) > 0) {
    gp1 <- gp1 + geom_curve(data = link_case_sub, curvature = -0.5, ncp = 10)
  }
  link_case_sub <- subset(link_case, stage == "11-14w" & sig & direction == 1)
  if(nrow(link_case_sub) > 0) {
    gp1 <- gp1 + geom_curve(data = link_case_sub, curvature = 0.5, ncp = 10)
  }
  gp1 <- gp1 + annotate(geom = "point", x = tss_pos, y = 0)
  # 19-22w
  gp2 <- gp0 + ylab("19-22w") #+ ggtitle(case_gene)
  link_case_sub <- subset(link_case, stage == "19-22w" & (! sig) & direction == -1)
  if(nrow(link_case_sub) > 0) {
    gp2 <- gp2 + geom_curve(data = link_case_sub, curvature = -0.5, ncp = 10, color = "grey90")
  }
  link_case_sub <- subset(link_case, stage == "19-22w" & (! sig) & direction == 1)
  if(nrow(link_case_sub) > 0) {
    gp2 <- gp2 + geom_curve(data = link_case_sub, curvature = 0.5, ncp = 10, color = "grey90")
  }
  link_case_sub <- subset(link_case, stage == "19-22w" & sig & direction == -1)
  if(nrow(link_case_sub) > 0) {
    gp2 <- gp2 + geom_curve(data = link_case_sub, curvature = -0.5, ncp = 10)
  }
  link_case_sub <- subset(link_case, stage == "19-22w" & sig & direction == 1)
  if(nrow(link_case_sub) > 0) {
    gp2 <- gp2 + geom_curve(data = link_case_sub, curvature = 0.5, ncp = 10)
  }
  gp2 <- gp2 + annotate(geom = "point", x = tss_pos, y = 0)
  # highlight dsl
  if(highlight.dsl) {
    gr <- geom_rect(data = dsl_case, aes(xmin = peakCenter - (xlim.max - xlim.min) * highlight.dsl.width, 
                                         xmax = peakCenter + (xlim.max - xlim.min) * highlight.dsl.width, 
                                         ymin = -Inf, ymax = Inf), fill = "purple", color = NA, alpha = 0.2)
    gp1 <- gp1 + gr
    gp2 <- gp2 + gr
  }
  # diff
  gpd <- ggplot(dsl_case) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") + 
    #geom_segment(aes(x = peakLeft, xend = peakRight, y = Correlation.x, yend = Correlation.x), color = "cornflowerblue", size = 1.5) + 
    #geom_segment(aes(x = peakLeft, xend = peakRight, y = Correlation.y, yend = Correlation.y), color = "limegreen", size = 1.5) + 
    geom_segment(data = dsl_case, aes(x = peakCenter, xend = peakCenter, y = Correlation.x, yend = Correlation.y, color = diff_abs > 0), 
                 arrow = arrow(length = unit(0.175, "inches")), show.legend = F) + 
    annotate(geom = "point", x = tss_pos, y = 0) + 
    scale_color_manual(values = c("blue", "red")) + 
    scale_x_continuous(limits = c(xlim.min, xlim.max)) + 
    scale_y_continuous(limits = c(min(dsl_case$Correlation.x, dsl_case$Correlation.y, -0.2), max(dsl_case$Correlation.x, dsl_case$Correlation.y)), 
                       breaks = seq(-1, 1, by = 0.2)) + 
    xlab(paste0("Genomic position", " (", unique(dsl_case$peakChr), ")")) + ylab("Correlation") + 
    theme(plot.margin = margin(7,margin.right,7,7, unit = "pt"))
  gp_LS <- list(gp1, gp2, gpd)
  return(gp_LS)
}
