
do_createDT <- function(expr_data, cellMetaData, do.norm = F, do.log = F, cell_num_cutoff = NULL) {
  if(! all(colnames(expr_data) == rownames(cellMetaData))) {
    stop("The cell sort in the two input are different, exit.")
  }
  if(do.norm) {
    cat("Normalize data...\n")
    expr_data <- sweep(expr_data, 2, colSums(expr_data), "/") * 1e6
  }
  if(do.log) {
    cat("Log transform data...\n")
    expr_data <- log2(expr_data + 1)
  }
  cellCluster_LS <- split(rownames(cellMetaData), cellMetaData$ident)
  # filtering
  if(! is.null(cell_num_cutoff)) {
    cellCluster_LS <- cellCluster_LS[lengths(cellCluster_LS) >= cell_num_cutoff]
  }
  
  res <- list()
  for(cell_type in names(cellCluster_LS)) {
    #cat(">", cell_type, "\n")
    res[[cell_type]][["data"]] <- expr_data[, colnames(expr_data) %in% cellCluster_LS[[cell_type]], drop = F]
  }
  return(res)
}

do_cor <- function(res_in, subsp = NULL, seed = 1, expr_cutoff = 0.1, mask = F, rm_outlier = F, method = "pearson", use = "all.obs") {
  cl <- makeCluster(min(length(res_in), 10), type = "FORK")
  expr_cor_LS <- parLapply(cl, res_in, function(x) {
    expr_input <- x[["data"]]
    # sub sampling
    if(! is.null(subsp)) {
      set.seed(seed)
      expr_input <- expr_input[, sample(1:ncol(expr_input), subsp)]
    }
    # filtering genes
    expr_input <- expr_input[rowSums(expr_input > 0) >= ncol(expr_input) * expr_cutoff, ]
    # mask zero
    if(mask) {
      expr_input[expr_input == 0] <- NA
    }
    # clac cor
    if(rm_outlier) {
      ngene <- nrow(expr_input)
      ncell <- ncol(expr_input)
      expr_cor <- matrix(NA, nrow = ngene, ncol = ngene)
      rownames(expr_cor) <- rownames(expr_input)
      colnames(expr_cor) <- rownames(expr_input)
      gene_mean <- rowMeans(expr_input)
      
      for(i in 1:ngene) {
        for(j in 1:ngene) {
          if(i==j) {
            expr_cor[i, j] <- 1
          } else {
            # identify outlier
            npro <- (expr_input[i, ] - gene_mean[i]) * (expr_input[j, ] - gene_mean[j])
            outlier <- which.max(npro)
            # recalculate mean without outlier
            if(method == "bicor") {
              expr_cor[i, j] <- WGCNA::cor(expr_input[i, -outlier], expr_input[j, -outlier], use = use)
            } else {
              expr_cor[i, j] <- cor(expr_input[i, -outlier], expr_input[j, -outlier], method = method, use = use)
            }
          }
        }
      }
    } else {
      if(method == "bicor") {
        expr_cor <- WGCNA::bicor(t(expr_input), use = use)
      } else {
        expr_cor <- cor(t(expr_input), method = method, use = use)
      }
    }
    return(expr_cor)
  })
  stopCluster(cl); rm(cl)
  for(cell_type in names(res_in)) {
    res_in[[cell_type]][["cor"]] <- expr_cor_LS[[cell_type]]
  }
  return(res_in)
}

do_findComBinding <- function(ctype, target_gene, tf_tf_aCt = 0.2, tf_tf_eCt = 0.2, tf_tg_eCt = 0.2, do.plot = F, do.write = F) {
  dir.create(path = paste0(OUT, "/res"), showWarnings = F, recursive = T)
  
  TF_and_peak <- read.table(file = paste0(OUT, "/TF_and_peak/", target_gene, ".txt"), header = F, sep = "\t", stringsAsFactors = F)
  colnames(TF_and_peak) <- c("tf", "peak")
  peak_id_sub <- rownames(peak_pos)[rownames(peak_pos) %in% unique(TF_and_peak$peak)]
  if(length(peak_id_sub) == 1) {
    cat("> Only one ATAC peak, skipped.\n")
    return(NULL)
  }
  
  ct <- ctype
  cell_type_label <- gsub(" ", "_", ct)
  if(! target_gene %in% rownames(res[[ct]]$cor)) {
    cat("> Target gene is not expressed in the specified cell type.\n")
    return(NULL)
  }
  
  acce_sub <- tes[[ct]]$data[peak_id_sub, ]
  cor_sub <- cor(t(acce_sub), method = "spearman")
  cor_sub[lower.tri(cor_sub, diag = T)] <- NA
  cor_sub_melted <- na.omit(melt(cor_sub))
  cor_sub_melted$Var1 <- as.character(cor_sub_melted$Var1)
  cor_sub_melted$Var2 <- as.character(cor_sub_melted$Var2)
  cor_sub_melted <- merge(cor_sub_melted, peak_pos, by.x = "Var1", by.y = 0)
  cor_sub_melted <- merge(cor_sub_melted, peak_pos, by.x = "Var2", by.y = 0)
  cor_sub_melted <- cor_sub_melted[, c(2,1,3:9)]
  #print(cor_sub_melted)
  
  # ggplot(cor_sub_melted, aes(Var1, Var2)) +
  #   geom_tile(aes(fill = value), color='white') +
  #   scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab') +
  #   theme(axis.text.x=element_text(angle=90), axis.ticks=element_blank(), axis.line=element_blank())
  
  if(do.plot) {
    x_label <- unique(c(cor_sub_melted$left.x, rev(cor_sub_melted$right.x)[1]))
    y_label <- unique(c(cor_sub_melted$left.y, rev(cor_sub_melted$right.y)[1]))
    tri_DF <- data.frame(x = c(min(x_label), min(x_label), max(y_label)), y = c(min(x_label), max(y_label), max(y_label)))
    gp <- ggplot() + 
      geom_polygon(data = tri_DF, aes(x = x, y = y), fill = "grey95") + 
      geom_rect(data = cor_sub_melted, mapping=aes(xmin = left.x, xmax = right.x, ymin = left.y, ymax = right.y, fill = value)) + 
      scale_fill_gradient2(low = 'blue', mid = "white", high = 'red', limits = c(-0.4, 0.4)) + 
      #scale_x_continuous(breaks = x_label, labels = x_label) + 
      #scale_y_continuous(breaks = y_label, labels = y_label) + 
      labs(fill = "Correlation") + 
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) + 
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
      theme(aspect.ratio = 1)
    
    ggsave(filename = paste0(OUT, "/res/", cell_type_label, "__", target_gene, ".pdf"), gp + theme(legend.position = "none"), width = 4, height = 4) 
    
    if(! file.exists(paste0(OUT, "/res/legend.pdf"))) {
      lg <- cowplot::get_legend(gp)
      ggsave(filename = paste0(OUT, "/res/legend.pdf"), lg, width = 1.25, height = 4)
    }
    
    cat("Region:", paste0(unique(cor_sub_melted$chr.x), ":", min(cor_sub_melted$left.x), "-", max(cor_sub_melted$right.y) + 14), "(-14)", "\n")
  }
  
  cor_sub_passed <- subset(cor_sub_melted, value > tf_tf_aCt)
  #print(cor_sub_passed)
  
  comBinding_LS <- apply(cor_sub_passed[, 1:3], 1, function(x) {
    tf_1 <- subset(TF_and_peak, peak == x[1], "tf", drop = T)
    tf_2 <- subset(TF_and_peak, peak == x[2], "tf", drop = T)
    tf_1 <- setdiff(tf_1, target_gene)
    tf_2 <- setdiff(tf_2, target_gene)
    cor_rna_sub <- res[[ct]]$cor[rownames(res[[ct]]$cor) %in% tf_1, colnames(res[[ct]]$cor) %in% tf_2, drop = F]
    cor_rna_sub_melted <- melt(cor_rna_sub)
    cor_rna_sub_melted$Var1 <- as.character(cor_rna_sub_melted$Var1)
    cor_rna_sub_melted$Var2 <- as.character(cor_rna_sub_melted$Var2)
    cor_rna_sub_melted <- subset(cor_rna_sub_melted, Var1 != Var2)
    cor_rna_sub_passed <- subset(cor_rna_sub_melted, value > tf_tf_eCt)
    cor_rna_sub_passed$target <- rep(target_gene, nrow(cor_rna_sub_passed))
    cor_rna_sub_passed$peak1 <- rep(x[1], nrow(cor_rna_sub_passed))
    cor_rna_sub_passed$peak2 <- rep(x[2], nrow(cor_rna_sub_passed))
    cor_rna_sub_passed$cor_acce <- rep(as.numeric(x[3]), nrow(cor_rna_sub_passed))
    cor_rna_sub_passed$cor_rna1 <- res[[ct]]$cor[cor_rna_sub_passed$Var1, target_gene]
    cor_rna_sub_passed$cor_rna2 <- res[[ct]]$cor[cor_rna_sub_passed$Var2, target_gene]
    cor_rna_sub_passed <- subset(cor_rna_sub_passed, cor_rna1 > tf_tg_eCt & cor_rna2 > tf_tg_eCt)
    colnames(cor_rna_sub_passed)[1:3] <- c("tf1", "tf2", "cor_rna0")
    return(cor_rna_sub_passed)
  })
  comBinding_DF <- do.call("rbind", comBinding_LS)
  comBinding_DF <- comBinding_DF[, c(1:4, 8:9, 5:7)]
  if(is.null(comBinding_DF)) {
    comBinding_DF <- as.data.frame(matrix(NA, nrow = 0, ncol = 9))
    colnames(comBinding_DF) <- c("tf1", "tf2", "cor_rna0", "target", "cor_rna1", "cor_rna2", "peak1", "peak2", "cor_acce")
  }
  if(do.write) {
    write.table(x = comBinding_DF, file = paste0(OUT, "/res/", cell_type_label, "__", target_gene, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  }
  
  # TFBS
  TFBS_in_peak <- read.table(file = paste0("03-expression/merged/comBinding/TFBS_in_peak/", target_gene, ".gtf"), header = F, sep = "\t", stringsAsFactors = F)
  TFBS_in_peak <- unique(TFBS_in_peak[, c(1,4,5,8,13)])
  colnames(TFBS_in_peak) <- c("chr", "left", "right", "tf", "peak")
  comBinding_tf_peak <- unique(data.frame(tf = c(comBinding_DF$tf1, comBinding_DF$tf2), 
                                          cor_rna = c(comBinding_DF$cor_rna1, comBinding_DF$cor_rna2), 
                                          peak = c(comBinding_DF$peak1, comBinding_DF$peak2), 
                                          stringsAsFactors = F))
  #print(head(comBinding_DF))
  #print(head(comBinding_tf_peak))
  comBinding_DF_withTFBS <- merge(comBinding_tf_peak, TFBS_in_peak, by = c("tf", "peak"), sort = F)
  
  gp <- ggplot(data = comBinding_DF_withTFBS) + geom_rect(aes(xmin = left, xmax = right, ymin = cor_rna - 0.005, ymax = cor_rna + 0.005)) + 
    ggrepel::geom_text_repel(aes(x = (left + right)/2, y = cor_rna - 0.015, label = tf), box.padding = 0.2, point.padding = 0) + 
    xlab(NULL) + ylab("TF-target correlation") + 
    scale_x_continuous(expand = c(0, 0), limits = c(min(cor_sub_melted$left.x), max(cor_sub_melted$right.y))) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if(do.plot) {
    ggsave(filename = paste0(OUT, "/res/", cell_type_label, "__", target_gene, "_TFBS.pdf"), gp, width = 6, height = 2.125)
  }
  return(comBinding_DF)
}

do_regression <- function(ctype, x_cb, return.model = F) {
  if(is.null(x_cb)) {
    return(NULL)
  }
  if(nrow(x_cb) == 0) {
    return(NULL)
  }
  
  ct <- ctype
  
  tf_id <- unique(c(x_cb$tf1, x_cb$tf2))
  target_id <- unique(x_cb$target)
  tf_other <- setdiff(intersect(TF_full_list, rownames(res[[ct]]$data)), c(tf_id, target_id))
  expr_sub <- as.data.frame(t(res[[ct]]$data[c(tf_id, tf_other, target_id), ]))

  # avoid specific char
  tf_id <- gsub("-", "___", tf_id)
  tf_other <- gsub("-", "___", tf_other)
  target_id <- gsub("-", "___", target_id)
  colnames(expr_sub) <- gsub("-", "___", colnames(expr_sub))
  
  # single TF
  y_s <- sapply(tf_id, function(x_tf) {
    reg_res <- lm(as.formula(paste(target_id, "~", x_tf)), data = expr_sub)
    y <- summary(reg_res)$adj.r.squared
  })
  
  # single TF with other TFs
  y_o <- sapply(tf_id, function(x_tf) {
    #cat(">", x_tf, "\n")
    res_rnd <- sapply(1:10, function(x_i) {
      #cat(">>", x_i, "\n")
      set.seed(x_i)
      tf_other_sub <- sample(tf_other, size = length(tf_id) - 1)
      #cat("  ", paste(target_id, "~", paste(c(x_tf, tf_other_sub), collapse = " + ")), "\n")
      reg_res <- lm(as.formula(paste(target_id, "~", paste(c(x_tf, tf_other_sub), collapse = " + "))), data = expr_sub)
      y <- summary(reg_res)$adj.r.squared
      return(y)
    })
    y <- res_rnd
  })
  y_o <- as.numeric(y_o)
  
  # all TFs
  reg_res <- lm(as.formula(paste(target_id, "~", paste(tf_id, collapse = " + "))), data = expr_sub)
  y_a <- summary(reg_res)$adj.r.squared
  pv <- pnorm(q = y_a, mean = mean(y_o), sd = sd(y_o), lower.tail = F)
  
  # plot(x = expr_sub[, target_id], y = fitted(reg_res), xlab = "Real expression", ylab = "Predicted expression", 
  #      xlim = c(min(expr_sub[, target_id], fitted(reg_res)), max(expr_sub[, target_id], fitted(reg_res))), 
  #      ylim = c(min(expr_sub[, target_id], fitted(reg_res)), max(expr_sub[, target_id], fitted(reg_res))))
  
  r2_adj_DF <- data.frame(target = target_id, variable = c(tf_id, "Single", "Single_withOther", "All"), 
                          value = c(y_s, mean(y_s), mean(y_o), y_a), p_value = pv, 
                          stringsAsFactors = F)
  rownames(r2_adj_DF) <- NULL
  
  # go back to specific char
  r2_adj_DF$target <- gsub("___", "-", r2_adj_DF$target)
  r2_adj_DF$variable <- gsub("___", "-", r2_adj_DF$variable)
  
  if(return.model) {
    return(reg_res)
  } else {
    return(r2_adj_DF)
  }
}
