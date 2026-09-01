# world infect sim

tbl_0 <- read.table('PATH TO LOG-TRANSFORMED RR TABLE', sep = '\t', check.names = FALSE)

probs_tbl_0 <- exp(tbl_0)
diag(probs_tbl_0) <- 0

all_regs <- colnames(tbl_0)
infected_regs <- rep(0, length(all_regs))

which(all_regs == 'Russia(EasternEurope)')
# which(all_regs == 'USA(NorthernAmerica)')

which_infected <- rep(0, length(all_regs))

inf_prob <- 0.008/(length(all_regs) - 1)
r0 <- 1.5

first_reg_array <- character()
inf_region_array <- matrix(nrow = length(all_regs) - 1, ncol = 5)
inf_prob_array <- matrix(nrow = length(all_regs) - 1, ncol = 5)
inf_times_array <- matrix(nrow = length(all_regs) - 1, ncol = 6)
inf_time_sd_array <- numeric()

for(first_reg in all_regs[all_regs != 'Russia(EasternEurope)']){
  j <- which(all_regs[all_regs != 'Russia(EasternEurope)'] == first_reg)  
  inf_region <- character()
    inf_time <- integer()
    for(i in 1:1000){
      for(interval_number in 0:200){
        if(interval_number == 0){
          infected_regs <- rep(0, length(all_regs))
          # infected_regs[which(all_regs == 'USA(NorthernAmerica)')] <- 1
          # infected_regs[which(all_regs == 'Turkey(WesternAsia)')] <- 1
          # infected_regs[which(all_regs == 'MiddleAfrica')] <- 1
          # infected_regs[which(all_regs == 'India(SouthernAsia)')] <- 1
          # infected_regs[which(all_regs == 'CentralAsia')] <- 1
          infected_regs[which(all_regs == first_reg)] <- 1
        }
        inf_tbl <- probs_tbl_0 * infected_regs * inf_prob
        inf_tbl_next <- apply(inf_tbl, 1, function(x){
          # print(x)
          as.numeric(runif(length(x)) < x)
          
        }) %>% t
        colnames(inf_tbl_next) <- rownames(inf_tbl_next)
        
        # REMOVE COMMENT HERE TO FORBID INFECTION FROM THE FIRST REGION TO RUSSIA
        # inf_tbl_next[which(all_regs == first_reg), 'Russia(EasternEurope)'] <- 0
        
        
        infected_regs <- infected_regs*r0 + colSums(inf_tbl_next)
        if(infected_regs[which(all_regs == 'Russia(EasternEurope)')] > 0){
          # print(all_regs[inf_tbl_next[, which(all_regs == 'Russia(EasternEurope)')] > 0])
          
          inf_region <- c(inf_region, all_regs[inf_tbl_next[, which(all_regs == 'Russia(EasternEurope)')] > 0])
          inf_time <- c(inf_time, interval_number)
          
          break
        }
        # print(interval_number)
        # print(infected_regs)
      }
      

    }
    print(first_reg)
    print((table(inf_region)/1000) %>% sort(TRUE))
    print(summary(inf_time))
    print(sd(inf_time))
    
    first_reg_array <- c(first_reg_array, first_reg)
    inf_region_array[j, ] <- names((table(inf_region)/1000) %>% sort(TRUE))[1:5]
    inf_prob_array[j, ] <- ((table(inf_region)/1000) %>% sort(TRUE))[1:5]
    inf_times_array[j, ] <- summary(inf_time)
    inf_time_sd_array <- c(inf_time_sd_array, sd(inf_time))
    
}

result_df <- cbind(first_reg_array, inf_region_array, inf_prob_array, inf_times_array, inf_time_sd_array)
colnames(result_df) <- c('first_reg_array', paste0('reg_', 1:5), paste0('reg_prob_', 1:5),
                         paste0('inf_time_', c('min', '1st.q', 'median', 'mean', '3rd.q', 'max', 'sd')))


write.table(result_df, 'PATH TO WRITE THE RESULTS', sep = '\t',
            row.names = FALSE, quote = FALSE)



