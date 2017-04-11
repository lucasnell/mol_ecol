# Make abundance files for different digestions

source('./wr_files/preamble.R')

rounded_p <- function(vec, round_digits = 4) {
    new_vec <- round(vec / sum(vec), round_digits)
    # If the rounded numbers don't sum to 1, then add the difference to the 
    # rounded number that differs most from the original, in the same direction as the
    # difference between sum(new_vec) and 1
    if (sum(new_vec) != 1) {
        newish_vec <- vec / sum(vec)
        diffs <- (newish_vec - new_vec) * 10^round_digits
        ind <- which(diffs == ifelse(sum(new_vec) < 1, min(diffs), max(diffs)))
        new_vec[ind] <- new_vec[ind] + (1 - sum(new_vec))
    }
    return(new_vec)
}

enz <- c('ApeKI', 'BstBI', 'BspEI')
n_samps <- 10
# Four scenarios:
scen_mat <- matrix(0, n_samps, 4)
scen_mat[,1] <- rep(0.1, 10)
scen_mat[,2] <- 1:10 %>% rounded_p
set.seed(656); scen_mat[,3] <- rlnorm(10) %>%  # <- More realistic given exp in growth terms
    rounded_p
scen_mat[,4] <- c(0.5, 0.25, 0.25, rep(0, 7))
# Turning to percents:
scen_mat <- scen_mat * 100

seq_names <- lapply(enz, 
                    function(e) {
                        f <- read_fasta(paste0('/Volumes/64gb/fasta/', e, '.fa.gz'))
                        names(f)
                    })


lapply(1:length(enz), 
       function(i) {
           
           .names <- seq_names[[i]]
           n_seqs_ps <- length(.names) / n_samps
           
           abund_mat <- lapply(1:ncol(scen_mat), 
                               function(i) rep(scen_mat[,i] / n_seqs_ps, 
                                               each = n_seqs_ps)) %>% 
               do.call(what = cbind, args = .)
           
           out_df <- as_data_frame(abund_mat) %>% 
               mutate_all(funs(format), scientific = FALSE) %>% 
               mutate(name = .names) %>% 
               select(name, everything())
           fn <- paste0('/Volumes/64gb/abundances/', enz[i], '.txt')
           write.table(out_df, fn, quote = FALSE, sep = ' ', row.names = FALSE, 
                       col.names = FALSE)
           
           return(invisible(NULL))
           
       })


