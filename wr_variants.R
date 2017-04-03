

# These are objects necessary for filtering sequences based on sequence size.

# If wr_preamble.R file hasn't been sourced, it needs to be.
# wr_preamble.R creates the .preamble_sourced object when it's run
if(!'.preamble_sourced' %in% ls(all.names = TRUE)) source('wr_preamble.R')



suppressPackageStartupMessages({
    library(gtools)
    library(parallel)
    library(Rcpp)
    library(RcppArmadillo)
})


sourceCpp('variants.cpp')


seg_sites <- 0.01414
seg_div <- 0.7072






pw_comp <- function(seq_vec) {
    seq_list <- .Internal(strsplit(seq_vec, '', FALSE, FALSE, FALSE))
    output <- sapply(seq_list, 
                     function(.s) {
                         pw_mat <- t(combn(.s, 2))
                         diffs <- ifelse(pw_mat[,1] == pw_mat[,2], 0, 1)
                         return(mean(diffs))
                     })
    return(output)
}




nt_freq <- function(N, divergence) {
    freq_mat <- combinations(N + 1, 4, 0:N, set = FALSE, repeats.allowed = TRUE)
    freq_sums <- rowSums(freq_mat)
    freq_mat <- freq_mat[freq_sums == N,]
    pw_divs <- pw_comp(apply(freq_mat, 1, function(x) cpp_make_seq(x)))
    min_diff_inds <- which(abs(pw_divs - divergence) == min(abs(pw_divs - divergence)))
    cat(sprintf('Minimum absolute difference from divergence = %s\n', 
                format(min(abs(pw_divs - divergence)), scientific = TRUE, digits = 4)))
    return(freq_mat[min_diff_inds,])
}






constr_objs <- function(dna_ss, seg_prop) {
    seqs <- as.character(dna_ss)
    seq_lens <- nchar(seqs)
    total_seg <- round(sum(seq_lens) * seg_prop, 0)
    rand_seqs <- sort(.Internal(sample(length(seq_lens), total_seg, replace = TRUE, 
                                       prob = seq_lens)))
    freq_len <- table(factor(rand_seqs, levels = 1:length(seq_lens))) %>% 
        as_data_frame %>% 
        mutate(len = seq_lens) %>% 
        select(n, len) %>% 
        as.matrix
    if (any(freq_len[,1] > freq_len[,2])) {
        stop('Retry with different seed.')
    }
    seq_obj <- list(`freq_len` = freq_len, `seqs` = seqs, `N` = length(seqs))
}





make_variants <- function(dna_ss, divergence, seg_prop, n_samps, cores = 1) {
    
    freq_mat <- nt_freq(n_samps, divergence)
    seq_obj <- constr_objs(dna_ss, seg_prop)
    
    .one <- function(i) {
        freq_len_row <- seq_obj$freq_len[i,]
        seq <- seq_obj$seqs[i]
        sites <- cpp_one_sites(freq_len_row)
        new_seqs <- cpp_change(seq, sites - 1, freq_mat, n_samps)
        return(new_seqs)
    }
    if (cores > 1) {
        new_sites <- mclapply(1:seq_obj$N, .one, mc.cores = cores)
    } else {
        new_sites <- lapply(1:seq_obj$N, .one)
    }
    new_sites <- c(new_sites, recursive = TRUE)
    dna_out <- DNAStringSet(new_sites)
    return(dna_out)
}


