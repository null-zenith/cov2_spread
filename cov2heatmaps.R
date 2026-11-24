make_log_rr_tbl <- function(tbl_path){
  tbl <- read.table(tbl_path, sep = '\t', header = TRUE, check.names = FALSE)
  
  tbl_log <- copy(tbl)

  for(i in 1:nrow(tbl)){
    for(j in 1:nrow(tbl)){
      if(i == j){
        tbl_log[i, j] <- 0
        next
        } else {
          value <- log(tbl[i, j]) + log(sum(tbl)) - log(sum(tbl[i, ])) - log(sum(tbl[, j]))
          tbl_log[i, j] <- ifelse(is.infinite(value), NA, value)
        }
    }
  }
  tbl_log
}



write.table(make_log_rr_tbl('/path/NO_BOOTSTRAP.tsv'),
            '/path/rr.tsv', sep = '\t', quote = FALSE)

library(pheatmap)
library(grid)
# library(ggplot2)
make_pheatmap <- function(rr_table_path, plot_name = '', nms = character(0), small.k = FALSE, continent.cols.path = NA, save.path = NA){
  rr_table <- read.table(rr_table_path, sep = '\t', header = TRUE, check.names = FALSE)
  
  if(length(nms) > 0){
    rr_table <- rr_table[nms, nms]
  }
  
  rr_table[] <- lapply(rr_table, function(x) {
    x[is.infinite(x)] <- NA
    x
  })
  
  row_dist <- as.dist(1 - cor(t(rr_table), use = "pairwise.complete.obs"))
  col_dist <- as.dist(1 - cor(rr_table, use = "pairwise.complete.obs"))
  

  if(!is.na(continent.cols.path)){
    col.table <- read.table(continent.cols.path, sep = '\t', header = TRUE, comment.char = '')
    rownames(col.table) <- col.table$Original
    
    p <- pheatmap(rr_table, fontsize_number = 20,
                  clustering_distance_rows = row_dist,
                  clustering_distance_cols = col_dist,
                  cutree_rows = ifelse(small.k, 3, 5), cutree_cols = ifelse(small.k, 3, 10), fontsize = 10, cellheight = 12, cellwidth = 12, main = plot_name, silent = TRUE)
    gt <- p$gtable
    
    row_order <- if (!is.null(p$tree_row)) p$tree_row$order else seq_len(nrow(rr_table))
    col_order <- if (!is.null(p$tree_col)) p$tree_col$order else seq_len(ncol(rr_table))
    
    # --- Recolor row names ---
    i <- which(gt$layout$name == "row_names")
    g <- gt$grobs[[i]]
    g$gp$col <- col.table[rownames(rr_table), "Continent_Color"][row_order]
    g$label    <- col.table[rownames(rr_table), "Short_Word_Form"][row_order]
    gt$grobs[[i]] <- g
    
    # --- Recolor column names ---
    j <- which(gt$layout$name == "col_names")
    g <- gt$grobs[[j]]
    g$gp$col <- col.table[rownames(rr_table), "Continent_Color"][row_order]
    g$label    <- col.table[rownames(rr_table), "Short_Word_Form"][row_order]
    gt$grobs[[j]] <- g
    
    
    
    if(!is.na(save.path)){
      ggsave(paste0(save.path, '/', gsub('tsv', 'png', basename(rr_table_path))), plot = gt, width = 15, height = 15, dpi = 300, bg = "white")
      
      
      m <- rr_table[row_order, col_order, drop = FALSE]
      rownames(m) <- gt$grobs[[ which(gt$layout$name == "row_names") ]]$label
      colnames(m) <- gt$grobs[[ which(gt$layout$name == "col_names") ]]$labe
      
      
      write.table(m, paste0(save.path, '/', basename(rr_table_path)), sep = '\t')
      
      
      
    } else {
      grid.newpage()
      grid.draw(gt)
    }
    
    
  } else {
    
    
    if(!is.na(save.path)){
      p <- pheatmap(rr_table, fontsize_number = 20,
                    clustering_distance_rows = row_dist,
                    clustering_distance_cols = col_dist,
                    cutree_rows = ifelse(small.k, 3, 5), cutree_cols = ifelse(small.k, 3, 10), fontsize = 10, cellheight = 12, cellwidth = 12, main = plot_name, silent = TRUE)
       ggsave(paste0(save.path, '/', gsub('tsv', 'png', basename(rr_table_path))), plot = p, width = 15, height = 15, dpi = 300, bg = "white")
      
      
      write.table(rr_table, paste0(save.path, '/', basename(rr_table_path)), sep = '\t')
      
      
    } else {
      pheatmap(rr_table, fontsize_number = 20,
               clustering_distance_rows = row_dist,
               clustering_distance_cols = col_dist,
               cutree_rows = ifelse(small.k, 3, 5), cutree_cols = ifelse(small.k, 3, 10), fontsize = 10, cellheight = 12, cellwidth = 12, main = plot_name)
    }
  }
}


make_pheatmap(rr_table_path = 'path',
              plot_name = 'Relative Risks for identical genomes',
              continent.cols.path = 'path',
              save.path = 'path')
