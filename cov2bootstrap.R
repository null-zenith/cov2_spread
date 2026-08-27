unss <- function(x, split = ',', fix = FALSE){
  splt <- strsplit(x, split = split)
  lst_lngth <- lapply(splt, length) %>% unlist
  if(fix){if(sum(lst_lngth == 0)) splt[lst_lngth == 0] <- ''}
  
  res <- unlist(splt)
  return(res)
}

suppressMessages(library(optparse))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
options(warn=1) 



option_list = list(
  make_option("--boot_tbl", action="store", default=NA, type='character',
              help="Path to the table for bootstrap"),
  make_option("--p1m1_tbl", action="store", default=NA, type='character',
              help="Path to p1m1 table, optional"),
  make_option("--N", action="store", default=1000, type='numeric',
              help="Number of permutations"),
  make_option("--Nsub", action="store", default=NA, type='numeric',
              help="Number of genomes to sample, default = NA (uses all)"),
  make_option("--out_folder", action="store", default='boot_tbl', type='character',
              help="Folder name of the bootstrap output, saved in the same forlder as boot table"))


opt = parse_args(OptionParser(option_list=option_list))



build_location_connections.v1 <- function(df) {
  
  df <- df %>%
    mutate(
      from.location = as.character(from.location),
      to.location = as.character(to.location)
    )
  
  all_ptrns <- unique(c(df$from, df$to))
  ptrn_lengths <- lengths(lapply(all_ptrns, unss))
  all_ptrns <- all_ptrns[order(ptrn_lengths)]
  all_locs <- sort(unique(c(df$from.location, df$to.location)))
  n_locs <- length(all_locs)
  
  mat <- matrix(0L, nrow = n_locs, ncol = n_locs,
                dimnames = list(all_locs, all_locs))
  
  loc_to_idx <- setNames(seq_along(all_locs), all_locs)
  
  df_from_split <- split(df$to.location, df$from)
  df_to_split <- split(df$from.location, df$to)
  
  df_ident <- df[df$from == df$to, ]
  ident_conn_by_ptrn <- split(df_ident$from.location, df_ident$from)
  
  print(paste0('Total number of ptrns: ', length(all_ptrns)))
  for(i in seq_along(all_ptrns)) {
    ptrn <- all_ptrns[i]
    
    if(i %% 10000 == 1) print(paste0('Analyzed ', i, ' ptrns...'))
    
    ident_conn <- ident_conn_by_ptrn[[ptrn]]
    if(is.null(ident_conn)) ident_conn <- character(0)
    
    from_locs <- df_from_split[[ptrn]]
    to_locs <- df_to_split[[ptrn]]
    
    all_conn <- unique(c(ident_conn, from_locs, to_locs))
    
    if(length(all_conn) < 2) next
    
    ident_idx <- loc_to_idx[ident_conn]
    all_idx <- loc_to_idx[all_conn]
    
    if(length(ident_idx) > 0) {
      mat[ident_idx, ident_idx] <- mat[ident_idx, ident_idx] + 1L
    }
    
    if(length(from_locs) > 0) {
      from_locs <- from_locs[from_locs != ptrn]
      if(length(from_locs) > 0) {
        from_idx <- loc_to_idx[from_locs]
        
        new_from_idx <- from_idx[!from_locs %in% ident_conn]
        
        if(length(new_from_idx) > 0) {
          mat[new_from_idx, ident_idx] <- mat[new_from_idx, ident_idx] + 1L
          mat[ident_idx, new_from_idx] <- mat[ident_idx, new_from_idx] + 1L
        }
      }
    }
    
    if(length(to_locs) > 0 && i > 1) {
      to_locs <- to_locs[to_locs != ptrn] # Remove self-connections
      if(length(to_locs) > 0) {
        valid_to <- !to_locs %in% all_ptrns[1:(i-1)]
        to_locs <- to_locs[valid_to]
        
        if(length(to_locs) > 0) {
          to_idx <- loc_to_idx[to_locs]
          
          new_to_idx <- to_idx[!to_locs %in% ident_conn]
          
          if(length(new_to_idx) > 0) {
            mat[new_to_idx, ident_idx] <- mat[new_to_idx, ident_idx] + 1L
            mat[ident_idx, new_to_idx] <- mat[ident_idx, new_to_idx] + 1L
          }
        }
      }
    }
  }
  
  diag(mat) <- 0L
  
  upper_tri <- upper.tri(mat, diag = TRUE)
  
  rows <- row(mat)[upper_tri]
  cols <- col(mat)[upper_tri]
  
  long <- data.frame(
    from = all_locs[rows],
    to   = all_locs[cols],
    n    = mat[upper_tri]
  )
  
  list(matrix = mat, long = long)
}





boot_tbl <- read.table(opt$boot_tbl, sep = '\t', header = TRUE, quote = '')

out_path <- ifelse(dirname(opt$boot_tbl) == '.',  paste0(opt$out_folder, '/') , paste0(dirname(opt$boot_tbl), '/', opt$out_folder, '/'))
print(paste0('Results are saved in ', out_path, ' folder'))
dir.create(out_path)

long_df <- data.frame(loc1 = character(), loc2 = character(), count = character(), iter = character())

from_to_tbl <- data.frame(from = character(0), to = character(0))

if(!is.na(opt$p1m1_tbl)){
  print('Found p1m1 table...')
  from_to_tbl <- read.table(opt$p1m1_tbl, sep = '\t', header = TRUE, quote = '')
  from_to_tbl <- from_to_tbl %>% separate_rows(to, sep = ";")
}



for(i in 1:opt$N){
  if ((i %% 10 == 1 & opt$N > 10) | (opt$N <= 10)) print(paste0('***************************************Starting iteration ', i))
  
  Nsub <- ifelse(is.na(opt$Nsub), nrow(boot_tbl), opt$Nsub)
  
  boot_rows <- sample(1:nrow(boot_tbl), Nsub, replace = TRUE)
  boot_tmp <- boot_tbl[boot_rows, ]
  
  boot_tmp <- boot_tmp[!is.na(boot_tmp$substitutions), ] %>% unique
  
  dup_subs <- unique(boot_tmp$substitutions[duplicated(boot_tmp$substitutions)])
  
  conn_tbl <- data.frame(from = dup_subs, to = dup_subs)
  
  
  if(!is.na(opt$p1m1_tbl)){
    conn_tbl <- rbind(conn_tbl, from_to_tbl)
    }
  

  setDT(conn_tbl)
  setDT(boot_tmp)
  

  setkey(conn_tbl, to)
  setkey(boot_tmp, substitutions)
  
  # First merge: inner join on 'to' (only keep rows where 'to' exists in boot_tmp$substitutions)
  mrgd <- unique(
    boot_tmp[conn_tbl, on = c(substitutions = "to"), nomatch = NULL, allow.cartesian = TRUE]
  )
  
  # Second merge: inner join on 'from' (only keep rows where 'from' exists in boot_tmp$substitutions)
  setkey(mrgd, from)
  mrgd <- unique(
    boot_tmp[mrgd, on = c(substitutions = "from"), nomatch = NULL, allow.cartesian = TRUE]
  )
  

  colnames(mrgd) <- c('from', 'from.location', 'to', 'to.location')
  mrgd <- mrgd[, c(1,3,2,4)]
  
  
  
  
  
  mrgd <- mrgd[mrgd$from.location != mrgd$to.location, ]
  
  rows_to_sort <- mrgd$from == mrgd$to
  
  
  # sorting regions for equal values of 'from' and 'to', because their order doesn't matter
  if (any(rows_to_sort)) {
    vals <- as.matrix(mrgd[rows_to_sort, c("from.location", "to.location")])
    sorted_vals <- t(apply(vals, 1, sort))
    sorted_vals <- as.data.frame(sorted_vals)
    # setDT(sorted_vals)
    mrgd[rows_to_sort, c("from.location", "to.location")] <- sorted_vals
  }
  mrgd <- unique(mrgd)
  
  
  all_levels <- sort(union(mrgd$from.location, mrgd$to.location))
  mrgd$from.location <- factor(mrgd$from.location, levels = all_levels)
  mrgd$to.location <- factor(mrgd$to.location, levels = all_levels)
  
  result <- build_location_connections.v1(mrgd)
  

  
  write.table(result$matrix, paste0(out_path, i, '.tsv'), sep = '\t', quote = FALSE)

  long_df <- rbind(long_df, data.frame(cbind(result$long, iteration = rep(i, nrow(result$long)))))
}
write.table(long_df, paste0(out_path, 'long_tbl.tsv'), sep = '\t', quote = FALSE)

print(paste0('Results are saved in ', out_path, ' folder'))

