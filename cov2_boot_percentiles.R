library(matrixStats)


dir_path <- ""



if (grepl("_01$", dir_path)) {
  base   <- sub("_01$", "", dir_path)
  dirs   <- paste0(base, sprintf("_%02d", 1:10))
  dirs   <- dirs[dir.exists(dirs)] 
  files  <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = '[0-9]+.tsv', full.names = TRUE)))
} else {
  files <- list.files(dir_path, pattern = '[0-9]+.tsv', full.names = TRUE)
}

read_one <- function(f) {
  as.matrix(read.table(f, row.names = 1, check.names = FALSE, sep = '\t'))
}

mats <- lapply(files, read_one)



# Build common label set (union across all replicates)

all_rn <- unique(unlist(lapply(mats, rownames)))
all_cn <- unique(unlist(lapply(mats, colnames)))

# symmetry requires the row and column label sets to match
if (!setequal(all_rn, all_cn)) {
  warning("Row and column label sets differ; using union of both.")
  all_rn <- all_cn <- union(all_rn, all_cn)
}
labels <- sort(all_rn)
L <- length(labels)


pad <- function(m) {
  full <- matrix(0, nrow = L, ncol = L, dimnames = list(labels, labels))
  full[rownames(m), colnames(m)] <- m
  full
}
mats <- lapply(mats, pad)
nr <- nc <- L; rn <- cn <- labels


# RR per replicate, with non-finite values -> NA
compute_RR <- function(m) {
  s      <- sum(m)
  rs     <- rowSums(m)
  cs     <- colSums(m)
  denom  <- outer(rs, cs)                  # rs[i] * cs[j]
  RR     <- (m * s) / denom                # 0/0 -> NaN, x/0 -> Inf
  RR[!is.finite(RR)] <- NA                 # handles zero row/col sums
  diag(RR) <- 1
  RR
}

RR_list <- lapply(mats, compute_RR)
RR_arr  <- simplify2array(RR_list)         # L x L x B
B <- dim(RR_arr)[3]
cat("stacked", B, "matrices into array", paste(dim(RR_arr), collapse = " x "), "\n")


# medians & quantiles 
flat     <- matrix(RR_arr, nrow = L * L, ncol = B)
med_flat <- rowMedians(flat, na.rm = TRUE)
qs_flat  <- rowQuantiles(flat,
                         probs = c(0.005, 0.025, 0.975, 0.995),
                         na.rm = TRUE)

# matrix
med_mat <- matrix(med_flat, nrow = L, ncol = L, dimnames = list(labels, labels))
diag(med_mat) <- 1
# saveRDS(med_mat, paste0(dir_path, "/RR_median_matrix.rds"))

# long tbl
idx <- which(upper.tri(med_mat, diag = FALSE), arr.ind = TRUE)
flat_idx <- idx[, 1] + (idx[, 2] - 1) * L   # row index into med_flat / qs_flat

long <- data.frame(
  name       = paste(labels[idx[, 1]], labels[idx[, 2]], sep = "|"),
  median     = med_flat[flat_idx],
  low.95.CI  = qs_flat[flat_idx, 2],
  high.95.CI = qs_flat[flat_idx, 3],
  low.99.CI  = qs_flat[flat_idx, 1],
  high.99.CI = qs_flat[flat_idx, 4],
  check.names = FALSE
)



sig99 <- (1 < long$low.99.CI | 1 > long$high.99.CI) &
  long$low.99.CI > 0 & long$high.99.CI > 0
sig95 <- (1 < long$low.95.CI | 1 > long$high.95.CI) &
  long$low.95.CI > 0 & long$high.95.CI > 0

long$signif <- ifelse(sig99, "**", ifelse(sig95, "*", ""))

write.table(long, paste0(dir_path, "/RR_summary_long.tsv"), row.names = FALSE, sep = '\t', quote = FALSE)

cat("done:", nrow(long), "cell pairs written\n")


# extra
make_symmetric_median_matrix <- function(long, sep = "|", diag_val = 1, log.output = TRUE) {
  # infer the full label set from the "row|col" names
  parts <- strsplit(long$name, sep, fixed = TRUE)
  labs  <- sort(unique(c(vapply(parts, `[`, "", 1),
                         vapply(parts, `[`, "", 2))))
  n <- length(labs)
  M <- matrix(NA_real_, n, n, dimnames = list(labs, labs))
  ri <- match(vapply(parts, `[`, "", 1), labs)
  ci <- match(vapply(parts, `[`, "", 2), labs)
  M[cbind(ri, ci)] <- long$median     # upper triangle only
  M[lower.tri(M)]  <- t(M)[lower.tri(M)]   # mirror to lower triangle
  diag(M) <- diag_val
  if(log.output){
    log(M)
    } else M
}

med_mat_from_long <- make_symmetric_median_matrix(long, log.output = TRUE)
write.table(med_mat_from_long, paste0(dir_path, "/RR_median_matrix_from_long_LOG.tsv"), sep = '\t', quote = FALSE)

med_mat_from_long <- make_symmetric_median_matrix(long, log.output = FALSE)
write.table(med_mat_from_long, paste0(dir_path, "/RR_median_matrix_from_long.tsv"), sep = '\t', quote = FALSE)

# saveRDS(med_mat_from_long, "RR_median_matrix_from_long.rds")


