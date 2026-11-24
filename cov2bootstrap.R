# unss <- function(x, split = ',', fix = FALSE){
#   splt <- strsplit(x, split = split)
#   lst_lngth <- lapply(splt, length) %>% unlist
#   if(fix){if(sum(lst_lngth == 0)) splt[lst_lngth == 0] <- ''}
#   
#   res <- unlist(splt)
#   return(res)
# }

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
              help="Number of stratified bootstrap replications, sampling with replacement; if N = 0, no bootstrap, N = 0 overrides Nsub"),
  make_option("--Nsub", action="store", default=NA, type='numeric',
              help="Number of genomes to sample, default = NA (uses all), if equal to numeric value, sampling is performed without replacement and stratification"),
  make_option("--out_folder", action="store", default='boot_tbl', type='character',
              help="Folder name of the bootstrap output, saved in the same forlder as boot table"))


opt = parse_args(OptionParser(option_list=option_list))



build_location_connections.v1 <- function(df) {
  dt <- as.data.table(df)
  dt[, `:=`(from.location = as.character(from.location),
            to.location   = as.character(to.location))]
  
  all_locs   <- sort(unique(c(dt$from.location, dt$to.location)))
  n_locs     <- length(all_locs)
  loc_to_idx <- setNames(seq_along(all_locs), all_locs)
  
  # --- identity locations, as integer indices ---------------------------
  ident <- unique(
    dt[from == to, .(from, loc_i = loc_to_idx[from.location])]
  )
  
  # --- anchor pairs: unordered pairs of identity locations per pattern ---
  # non-equi on INTEGER columns (character not supported)
  anchor_pairs <- if (nrow(ident) > 0L) {
    ident[, .(from, a = loc_i)][
      ident[, .(from, b = loc_i)],
      on              = .(from, a < b),
      j               = .(from = x.from, a = x.a, b = i.b),
      nomatch         = NULL,
      allow.cartesian = TRUE
    ][, .(a, b)]
  } else {
    data.table(a = integer(), b = integer())
  }
  
  # --- destinations as integer indices ----------------------------------
  dest <- unique(dt[, .(from, loc_i = loc_to_idx[to.location])])
  # drop destinations that are already identity locations for the same pattern
  dest <- dest[!ident, on = .(from, loc_i)]
  
  # --- star: each destination × each identity location per pattern ------
  # equality on `from` only -> character key on `from` is fine here
  star_pairs <- if (nrow(dest) > 0L && nrow(ident) > 0L) {
    merge(
      dest[, .(from, a = loc_i)],
      ident[, .(from, b = loc_i)],
      by              = "from",
      allow.cartesian = TRUE
    )[a != b, .(a, b)]
  } else {
    data.table(a = integer(), b = integer())
  }
  
  # --- combine, canonicalize, count -------------------------------------
  pairs <- rbindlist(list(anchor_pairs, star_pairs))
  pairs <- pairs[a != b]
  
  if (nrow(pairs) > 0L) {
    pairs[, `:=`(lo = pmin(a, b), hi = pmax(a, b))]
    counts <- pairs[, .N, by = .(lo, hi)]
  } else {
    counts <- data.table(lo = integer(), hi = integer(), N = integer())
  }
  
  # --- scatter symmetrically -------------------------------------------
  mat <- matrix(0L, n_locs, n_locs, dimnames = list(all_locs, all_locs))
  if (nrow(counts) > 0L) {
    mat[cbind(counts$lo, counts$hi)] <- counts$N
    mat[cbind(counts$hi, counts$lo)] <- counts$N
  }
  
  upper_tri <- upper.tri(mat, diag = TRUE)
  rows <- row(mat)[upper_tri]
  cols <- col(mat)[upper_tri]
  long <- data.frame(from = all_locs[rows], to = all_locs[cols],
                     n    = mat[upper_tri])
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


if(opt$N == 0){
  print('No bootstrap mode running...')
  
  boot_tmp <- copy(boot_tbl)
  
  boot_tmp <- boot_tmp[!is.na(boot_tmp$substitutions), ] %>% unique
  

  conn_tbl <- data.frame(from = copy(boot_tmp$substitutions), to = copy(boot_tmp$substitutions))
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
  
  
  rows_to_sort <- mrgd$from == mrgd$to
  
  
  # sorting regions for equal values of 'from' and 'to', because their order doesn't matter
  if (any(rows_to_sort)) {
    vals <- as.matrix(mrgd[rows_to_sort, c("from.location", "to.location")])
    sorted_vals <- t(apply(vals, 1, sort))
    sorted_vals <- as.data.frame(sorted_vals)
    mrgd[rows_to_sort, c("from.location", "to.location")] <- sorted_vals
  }
  mrgd <- unique(mrgd)
  
  
  all_levels <- sort(union(mrgd$from.location, mrgd$to.location))
  mrgd$from.location <- factor(mrgd$from.location, levels = all_levels)
  mrgd$to.location <- factor(mrgd$to.location, levels = all_levels)
  
  result <- build_location_connections.v1(mrgd)
  
  
  write.table(result$matrix, paste0(out_path, 'NO_BOOTSTRAP.tsv'), sep = '\t', quote = FALSE)
  # count_mtrx <- as.data.frame(count_mtrx)
  long_df <- data.frame(cbind(result$long))
  
} else {
  for(j in 1:opt$N){
    if ((j %% 10 == 1 & opt$N > 10) | (opt$N <= 10)) print(paste0('***************************************Starting iteration ', j))
    
    if(is.na(opt$Nsub)){
      print('Stratified sampling mode running...')
      boot_tmp <- boot_tbl %>%
        group_by(across(2)) %>%
        slice_sample(prop = 1, replace = TRUE) %>%
        ungroup()
    } else {
      print('Subsampling mode running...')
      boot_tmp <- boot_tbl %>%
        slice_sample(n = opt$Nsub, replace = FALSE) %>%
        ungroup()
      
      # boot_rows <- sample(1:nrow(boot_tbl), opt$Nsub, replace = TRUE)
      # print(paste0('boot rows:', Nsub))
      # boot_tmp <- boot_tbl[boot_rows, ]
      
    }

    boot_tmp <- boot_tmp[!is.na(boot_tmp$substitutions), ] %>% unique

    conn_tbl <- data.frame(from = copy(boot_tmp$substitutions), to = copy(boot_tmp$substitutions))
    
    
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
    
    
    rows_to_sort <- mrgd$from == mrgd$to
    
    
    # sorting regions for equal values of 'from' and 'to', because their order doesn't matter
    if (any(rows_to_sort)) {
      vals <- as.matrix(mrgd[rows_to_sort, c("from.location", "to.location")])
      sorted_vals <- t(apply(vals, 1, sort))
      sorted_vals <- as.data.frame(sorted_vals)
      mrgd[rows_to_sort, c("from.location", "to.location")] <- sorted_vals
    }
    mrgd <- unique(mrgd)
    
    
    all_levels <- sort(union(mrgd$from.location, mrgd$to.location))
    mrgd$from.location <- factor(mrgd$from.location, levels = all_levels)
    mrgd$to.location <- factor(mrgd$to.location, levels = all_levels)
    
    result <- build_location_connections.v1(mrgd)
    
    
    write.table(result$matrix, paste0(out_path, j, '.tsv'), sep = '\t', quote = FALSE)
    # count_mtrx <- as.data.frame(count_mtrx)
    long_df <- rbind(long_df, data.frame(cbind(result$long, iteration = rep(j, nrow(result$long)))))
  }
}

if(opt$N == 0){
  write.table(long_df, paste0(out_path, 'NO_BOOTSTRAP_long_tbl.tsv'), sep = '\t', quote = FALSE)
} else write.table(long_df, paste0(out_path, 'long_tbl.tsv'), sep = '\t', quote = FALSE)

print(paste0('Results are saved in ', out_path, ' folder'))

