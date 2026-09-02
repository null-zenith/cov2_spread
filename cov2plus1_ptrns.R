suppressMessages(library(magrittr))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
options(warn=1)



option_list = list(
  make_option("--nc_path", action="store", default=NA, type='character',
              help="Path to the Nextclade table(s), comma-separated"),
  make_option("--skip", action="store", default=FALSE, type='logical',
              help="If TRUE, skip to plus1 part"))

opt = parse_args(OptionParser(option_list=option_list))

unss <- function(x, split = ',', fix = FALSE){
  splt <- strsplit(x, split = split)
  lst_lngth <- lapply(splt, length) %>% unlist
  if(fix){if(sum(lst_lngth == 0)) splt[lst_lngth == 0] <- ''}
  
  res <- unlist(splt)
  return(res)
}


if(!opt$skip){

  print('Filtering...')
  nc_tbl <- read.table(opt$nc_path, sep = '\t', header = TRUE, quote = '', comment.char = '')
  nc_tbl <- nc_tbl[(nc_tbl$totalMissing <= 300 & nc_tbl$qc.overallStatus == 'good' & (nc_tbl$alignmentEnd - nc_tbl$alignmentStart > 29000)), c('seqName', 'Nextclade_pango', 'partiallyAliased', 'substitutions')]
  write.table(nc_tbl, gsub('tsv', 'Nfiltered.tsv', opt$nc_path), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

  nc_tbl <- nc_tbl[, c('Nextclade_pango', 'partiallyAliased', 'substitutions')]
  
  print('Finding discordant subs and lineages... Discordant records are reported only and NOT filtered out.')
  nc_tbl_sub_diff_pango <-  nc_tbl %>% group_by(substitutions) %>% filter(n_distinct(partiallyAliased) > 1) %>% distinct(substitutions, partiallyAliased) %>% ungroup()
  nc_tbl_sub_diff_pango <- nc_tbl_sub_diff_pango[order(nc_tbl_sub_diff_pango$partiallyAliased), ]
  write.table(nc_tbl_sub_diff_pango[, c("partiallyAliased", "substitutions")], gsub('tsv', 'diff_pango.tsv', opt$nc_path), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  
  
  nc_tbl <- unique(nc_tbl)

  # Removing duplicates
  nc_tbl$seqName_short <- gsub('\\|.*', '', nc_tbl$seqName)
  dup_names <- nc_tbl$seqName_short[duplicated(nc_tbl$seqName_short)] %>% unique
  nc_tbl <- nc_tbl[!(nc_tbl$seqName_short %in% dup_names), ]

  write.table(nc_tbl, gsub('tsv', 'short.tsv', opt$nc_path), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  nsub <- sapply(nc_tbl$substitutions, function(x) length(unss(x)))
  nsub_list <- sort(unique(nsub))
  print(paste0('number of mutations up to ', max(nsub)))
  
  system(paste0('mkdir ', dirname(opt$nc_path), '/by_nsub/'))
  for(i in nsub_list){
    print(i)
    write.table(nc_tbl[nsub == i, ], paste0(dirname(opt$nc_path), '/by_nsub/', i, '.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
  }
}



print('Finding plus one ...')
CON <- file(paste0(dirname(opt$nc_path), '/plus1.txt'), 'a')
print(CON)
fls <- system(paste0('ls ', dirname(opt$nc_path), '/by_nsub/*.tsv'), intern = TRUE)

nsub <- as.numeric(gsub('.tsv', '', basename(fls)))

for(i in 0:max(nsub)){
  print(paste0('Processing file with ', i, ' mutations...'))
  # print(paste0(dirname(opt$nc_path), '/by_nsub/', i, '.tsv'))
  if (file.exists(paste0(dirname(opt$nc_path), '/by_nsub/', i, '.tsv')) & file.exists(paste0(dirname(opt$nc_path), '/by_nsub/', i+1, '.tsv'))){
    tbl1 <- read.table(paste0(dirname(opt$nc_path), '/by_nsub/', i, '.tsv'), sep = '\t', header = TRUE)
    tbl1[is.na(tbl1)] <- ''
    uniq_ptrns1 <- tbl1$substitutions %>% unique

    tbl2 <- read.table(paste0(dirname(opt$nc_path), '/by_nsub/', i+1, '.tsv'), sep = '\t', header = TRUE)
    uniq_ptrns2 <- tbl2$substitutions %>% unique

    if(i == 0){
      write(paste0('root -> ', paste0(uniq_ptrns2, collapse = ';')), CON)
    } else {
      uniq_subs2 <- unss(uniq_ptrns2, split = ',') %>% unique
      rev_index <- sapply(uniq_subs2, function(x) uniq_ptrns2[grepl(x, uniq_ptrns2)]) %>% as.list

      # print(uniq_ptrns1)
      for(j in 1:length(uniq_ptrns1)){

        current_muts <- unss(uniq_ptrns1[j])

        for(k in 1:length(current_muts)){
          if (k == 1){
            ptrns_tmp <- rev_index[[current_muts[k]]]
          } else {
            ptrns_tmp <- ptrns_tmp[ptrns_tmp %in% rev_index[[current_muts[k]]]]
            if(length(ptrns_tmp) == 0) break
          }


        }
        if(length(ptrns_tmp) != 0) write(paste0(uniq_ptrns1[j], ' -> ', paste0(ptrns_tmp, collapse = ';')), CON) # print(paste0(uniq_ptrns1[j], ' -> ', paste0(ptrns_tmp, collapse = ';'))) 

      }
      # print(ptrns_tmp)

    }
  }
}
close(CON)

                          
p1m1_world <- readLines(CON)

p1m1_df <- data.frame(from = character(length(p1m1_world)), to = character(length(p1m1_world)))

p1m1_df$from <- gsub(' -> .*', '', p1m1_world)
p1m1_df$to <- gsub('.* -> ', '', p1m1_world)


##################################################################################################FILTRATION OF AMBIGUOUS p1m1 CONNECTIONS
p1m1_df_to <- unss(p1m1_df$to, split = ';')
p1m1_df_to_dup <- p1m1_df_to[duplicated(p1m1_df_to)] %>% unique

p1m1_df$to <- sapply(p1m1_df$to, function(x){
  splt <- unss(x, split = ';')
  paste0(splt[!(splt %in% p1m1_df_to_dup)], collapse = ';')
})

# dropping p1m1 pairs where no "to" genome left
p1m1_df <- p1m1_df[p1m1_df$to != '', ]

write.table(p1m1_df, gsub('.txt', '_from_to.tsv', CON), sep = '\t', row.names = FALSE, quote = FALSE)
