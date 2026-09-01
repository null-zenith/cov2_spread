library(data.table)
library(magrittr)
library(pheatmap)
library(dplyr)
library(ggplot2)

aggr_boot_conn <- function(table_path, pseudocount = 0, remove.zero = TRUE){


  tbl <- read.table(table_path, header = TRUE, sep = '\t')
  
  
  ######################################### AGGREGATIION BLOCK ################################################
  # if you decide to join some of the regions, then it's easier to do it now than to make a new boot table and run bootsrapping again
  # aggregation ued in paper is further below
  
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

#########################################################################################    END OF AGGRERATION






  setDT(tbl)
  
  locations <- unique(c(tbl$from, tbl$to)) %>% sort
  
  tbl$n[tbl$n == 0] <- pseudocount
  

  tbl[, total_n_iter := sum(n), by = iteration]
  nodes_dt <- unique(tbl[, .(node = c(from, to)), by = iteration])

  node_sums <- nodes_dt[, .(sum_n = {
    tbl[iteration == .BY$iteration & (from == node | to == node), sum(n)]
  }), by = .(iteration, node)]
  tbl[node_sums, sum_n_from := i.sum_n, on = .(iteration, from = node)]
  tbl[node_sums, sum_n_to   := i.sum_n, on = .(iteration, to   = node)]

  

  tbl[, rr := ifelse(from == to, 0, 
                     2 * (as.numeric(tbl[[3]]) + pseudocount) * total_n_iter / (sum_n_from * sum_n_to))]


  nodes <- sort(unique(c(tbl$from, tbl$to)))
  n_nodes <- length(nodes)
  total_iter <- uniqueN(tbl$iteration)
  tbl[, `:=`(node1 = pmin(from, to), node2 = pmax(from, to))]
  pairs_dt <- tbl[node1 != node2, .(
    rr_values = list(rr),
    n_iter    = uniqueN(iteration)
  ), by = .(node1, node2)]
  
  ###################################CALCULATION OF STAR MATRIX: Checking if the range of bootstrapped RR values crosses 1; if not, this cell gets a star
  pairs_dt[, value := {
    if (n_iter < total_iter / 2) {
      NA_real_
    } else {
      rr_vec <- unlist(rr_values)
      if (min(rr_vec) <= 1 && max(rr_vec) >= 1) 1 else mean(rr_vec, na.rm = TRUE)
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
  
  write.table(star_matrix, gsub('.tsv', '.starmatrix.tsv', table_path), sep = '\t', quote = FALSE)
  ###################################END OF STAR MATRIX BLOCK
  
  
  ###################################CALCULATION OF MEAN RR table
  pairs_dt <- tbl[node1 != node2, .(
    rr_values = list(rr),
    n_iter    = uniqueN(iteration)
  ), by = .(node1, node2)]
  
  pairs_dt[, value := {
    if (n_iter < total_iter / 2) {
      NA_real_
    } else {
      rr_vec <- unlist(rr_values)
      mean(rr_vec, na.rm = TRUE)
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

  
  
  # filter out obviously empty rows and columns
  rows_to_keep <- !apply(mat, 1, function(row) all(is.na(row) | row == 0 | row == 1))
  all(rows_to_keep)
  mat_sym <- mat[rows_to_keep, rows_to_keep, drop = FALSE]
  mat_sym_log <- log(mat_sym)
  
  diag(mat_sym_log) <- 0
  mat_sym_log[is.infinite(mat_sym_log)] <- NA
  
  write.table(mat_sym_log, gsub('.tsv', '.RR.tsv', table_path), sep = '\t', quote = FALSE)
  ##################################END OF RR CALCULATION
  
  # quick visualization
  row_dist <- as.dist(1 - cor(t(mat_sym_log), use = "pairwise.complete.obs"))
  col_dist <- as.dist(1 - cor(mat_sym_log, use = "pairwise.complete.obs"))
  
  pheatmap(mat_sym_log, display_numbers = star_matrix, fontsize_number = 18, 
           clustering_distance_rows = row_dist,
           clustering_distance_cols = col_dist,
           cutree_rows = 5, cutree_cols = 10, fontsize = 8, cellheight = 12, cellwidth = 12, main = "Corr dist, Russia, p1m1, 1000")

}




# example
# aggr_boot_conn('PATH')