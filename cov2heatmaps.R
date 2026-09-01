library(pheatmap)

make_star_matrix <- function(table_path){
  
  tbl <- read.table(table_path, header = TRUE, sep = '\t')
  
  
  #################################### RUS SPLIT ONLY
  
  rus_aggregation <- c(
    "SouthernFD(Russia)" = "Southern+NorthCaucasianFD(Russia)",
    "NorthCaucasianFD(Russia)" = "Southern+NorthCaucasianFD(Russia)")
  
  tbl <- tbl %>%
    group_by(iteration) %>%
    mutate(
      from = ifelse(from %in% names(rus_aggregation), rus_aggregation[from], from),
      to = ifelse(to %in% names(rus_aggregation), rus_aggregation[to], to)
    ) %>%
    ungroup() %>%
    mutate(
      # Sort from and to alphabetically within each row
      from_new = pmin(from, to),
      to_new   = pmax(from, to),
      from = from_new,
      to   = to_new
    ) %>%
    select(-from_new, -to_new) %>%
    group_by(iteration, from, to) %>%
    summarise(n = sum(n), .groups = "drop")
  
  tbl <- tbl[, c(2:4, 1)]
  colnames(tbl) <- c('from', 'to', 'n', 'iteration')
  
  
  ######################################################################################### WORLD
  
  PI_aggregation <- c(
    "Polynesia" = "PacificIslands",
    "Melanesia" = "PacificIslands",
    "Micronesia" = "PacificIslands")
  
  tbl <- tbl %>%
    group_by(iteration) %>%
    mutate(
      from = ifelse(from %in% names(PI_aggregation), PI_aggregation[from], from),
      to = ifelse(to %in% names(PI_aggregation), PI_aggregation[to], to)
    ) %>%
    ungroup() %>%
    mutate(
      # Sort from and to alphabetically within each row
      from_new = pmin(from, to),
      to_new   = pmax(from, to),
      from = from_new,
      to   = to_new
    ) %>%
    select(-from_new, -to_new) %>%
    group_by(iteration, from, to) %>%
    summarise(n = sum(n), .groups = "drop")
  
  tbl <- tbl[, c(2:4, 1)]
  colnames(tbl) <- c('from', 'to', 'n', 'iteration')
  
  ######################################################################################### SEA ONLY
  
  SEA_aggregation <- c(
    "Brunei(SouthEasternAsia)" = "SouthEasternAsia"
  )
  
  tbl <- tbl %>%
    group_by(iteration) %>%
    mutate(
      from = ifelse(from %in% names(SEA_aggregation), SEA_aggregation[from], from),
      to = ifelse(to %in% names(SEA_aggregation), SEA_aggregation[to], to)
    ) %>%
    ungroup() %>%
    mutate(
      # Sort from and to alphabetically within each row
      from_new = pmin(from, to),
      to_new   = pmax(from, to),
      from = from_new,
      to   = to_new
    ) %>%
    select(-from_new, -to_new) %>%
    group_by(iteration, from, to) %>%
    summarise(n = sum(n), .groups = "drop")
  
  tbl <- tbl[, c(2:4, 1)]
  colnames(tbl) <- c('from', 'to', 'n', 'iteration')
  
  #########################################################################################  
  
  setDT(tbl)
  
  # locations <- unique(c(tbl$from, tbl$to)) %>% sort
  tbl[, total_n_iter := sum(n), by = iteration]
  nodes_dt <- unique(tbl[, .(node = c(from, to)), by = iteration])
  node_sums <- nodes_dt[, .(sum_n = {
    tbl[iteration == .BY$iteration & (from == node | to == node), sum(n)]
  }), by = .(iteration, node)]
  tbl[node_sums, sum_n_from := i.sum_n, on = .(iteration, from = node)]
  tbl[node_sums, sum_n_to   := i.sum_n, on = .(iteration, to   = node)]
  tbl[, rr := ifelse(from == to, 0, 
                     2 * as.numeric(n) * as.numeric(total_n_iter) / (as.numeric(sum_n_from) * as.numeric(sum_n_to)))]
  
  nodes <- sort(unique(c(tbl$from, tbl$to)))
  n_nodes <- length(nodes)
  total_iter <- uniqueN(tbl$iteration)
  tbl[, `:=`(node1 = pmin(from, to), node2 = pmax(from, to))]
  pairs_dt <- tbl[node1 != node2, .(
    rr_values = list(rr),                # list of all rr for this pair
    n_iter    = uniqueN(iteration)       # distinct iterations where pair appears
  ), by = .(node1, node2)]
  
  
  
  # for(i in 1:nrow(pairs_dt)){
  #   if(pairs_dt$n_iter[i] < total_iter / 2){
  #     pairs_dt$value[i] <- 1
  #   } else {
  #     rr_vec <- unlist(pairs_dt$rr_values[i])
  #     pairs_dt$value[i] <- ifelse((min(rr_vec, na.rm = TRUE) <= 1 && max(rr_vec, na.rm = TRUE) >= 1) | (sum(rr_vec == 0, na.rm = TRUE) + sum(is.nan(rr_vec) > 0.5*total_iter)), 1, mean(rr_vec, na.rm = TRUE))
  #     }
  # }
  
  
  pairs_dt[, value := {
    if (n_iter < total_iter / 2) {
      1
    } else {
      rr_vec <- unlist(rr_values)
      # mean(rr_vec, na.rm = TRUE)
      ############ median(rr_vec, na.rm = TRUE)
      # q <- quantile(rr_vec, probs = c(0.025, 0.975), na.rm = TRUE)
      # if (q[1] <= 1 && q[2] >= 1) 1 else mean(rr_vec, na.rm = TRUE)
      
      
      # if (min(rr_vec) <= 1 && max(rr_vec) >= 1) 1 else mean(rr_vec, na.rm = TRUE)
      if ((min(rr_vec, na.rm = TRUE) <= 1 && max(rr_vec, na.rm = TRUE) >= 1) | (sum(rr_vec == 0, na.rm = TRUE) + sum(is.nan(rr_vec)) > 0.5*total_iter)) 1 else mean(rr_vec, na.rm = TRUE)
    }
  }, by = .(node1, node2)]
  
  mat <- matrix(NA, nrow = n_nodes, ncol = n_nodes,
                dimnames = list(nodes, nodes))
  diag(mat) <- 1
  
  
  for (i in 1:nrow(pairs_dt)) {
    r <- as.character(pairs_dt$node1[i])
    c <- as.character(pairs_dt$node2[i])
    val <- pairs_dt$value[i]
    mat[r, c] <- val
    mat[c, r] <- val
  }
  
  
  star_matrix <- ifelse((mat != 1), "\u2217", "")
  
  write.table(star_matrix, paste0(table_path, '.CI'), sep = '\t', quote = FALSE)
  
}

plot_heatmap <- function(rr_path, star_path, add_spaces = FALSE, save.tbl = FALSE, plot_name = "", drop_names = FALSE){
  
  rr_table <- read.table(rr_path, header = TRUE, sep = '\t', check.names = FALSE)
  star_table <- read.table(star_path, header = TRUE, sep = '\t', check.names = FALSE, colClasses = "character")
  
  add.spaces <- function(string){
    if(grepl('^USA', string)) return(string)
    if(grepl('Australiaand', string)) string <- gsub('Australiaand', 'Australia and', string)
    string1 <- gsub("(?<![(])([A-Z(])", " \\1", string, perl = TRUE)
    sub("^ ", "", string1)
  }
  
  if(drop_names){
    
    star_table <- star_table[rownames(rr_table), colnames(rr_table)]
    
  }
  
  if(any(rownames(rr_table) != colnames(rr_table))) stop('Error in rr_table names')
  if(any(rownames(star_table) != colnames(star_table))) stop('Error in star_table names')
  if(any(rownames(rr_table) != rownames(star_table))) stop('Inconsistency between rr_table and star_table names')
  
  
  save_dir <- dirname(rr_path) 
  if(add_spaces){
    save_dir <- paste0(dirname(rr_path), '/FORMATTED')
    dir.create(save_dir)
    for(i in 1:length(colnames(rr_table))) colnames(rr_table)[i] <- add.spaces(colnames(rr_table)[i])
    for(i in 1:length(rownames(rr_table))) rownames(rr_table)[i] <- add.spaces(rownames(rr_table)[i])
    
    
    
    for(i in 1:length(colnames(star_table))) colnames(star_table)[i] <- add.spaces(colnames(star_table)[i])
    for(i in 1:length(rownames(star_table))) rownames(star_table)[i] <- add.spaces(rownames(star_table)[i])
  }
  print(save_dir)
  
  if(save.tbl){
    write.table(rr_table, paste0(save_dir, '/', basename(rr_path)), sep = '\t', quote = FALSE)
    write.table(star_table, paste0(save_dir, '/', basename(star_path)), sep = '\t', quote = FALSE)
  }
  
  
  row_dist <- as.dist(1 - cor(t(rr_table), use = "pairwise.complete.obs"))
  col_dist <- as.dist(1 - cor(rr_table, use = "pairwise.complete.obs"))
  phmp <- pheatmap(rr_table, display_numbers = star_table, fontsize_number = 10,
                   clustering_distance_rows = row_dist,
                   clustering_distance_cols = col_dist,
                   cutree_rows = 5, cutree_cols = 10, fontsize = 8, cellheight = 12, cellwidth = 12, main = plot_name)
  
  ggsave(paste0(save_dir, '/', gsub('tsv', 'png', basename(rr_path))), plot = phmp, width = 15, height = 15, dpi = 300, bg = "white")
  
  
}





make_star_matrix(table_path = '...')



plot_heatmap(rr_path = '...',
             star_path = '...',
             plot_name = '...')


