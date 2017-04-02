#' ---
#' title: 'Create structural variants from reference genome'
#' author: 'Lucas Nell'
#' date: '24 March 2017'
#' output: 
#'   pdf_document:
#'     toc: true
#'     number_sections: true
#' ---
#' 
#' *Updated `r format(Sys.Date(), '%d %B %Y')`*
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages({
    library(magrittr)
    library(ggplot2)
    library(purrr)
    library(dplyr)
    library(ShortRead)
    library(gtools)
    library(parallel)
    library(Rcpp)
    library(RcppArmadillo)
})
#' 
#' Some of this code is written in C++, so I need to load that file.
#' 
#+ load_variants
sourceCpp('variants.cpp')
#' 
#' 
#' 
#+ set_theme, echo = FALSE
# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))
#' 
#' 
#' 
#' 
#' 
#' 
#' # Data on which to base simulations
#' 
#' Reference:
#' 
#' \hangindent=0.25in Bickel, R. D., J. P. Dunham, and J. A. Brisson. 2013. 
#' Widespread selection across coding and noncoding DNA in the pea aphid genome. 
#' *G3: Genes|Genomes|Genetics* __3__:993–1001.
#' Available from 
#' [http://www.g3journal.org/content/3/6/993](http://www.g3journal.org/content/3/6/993)
#' 
#' 
#' The main points are below (all quotes are from p 996):
#' 
#' * "We sequenced 21 genetically distinct lines of pea aphids..."
#' * "...we calculated $F_{st}$ levels across the genome, comparing 11
#'   pea aphid lines from the Northeast US (New York and Massachusetts)
#'   and 10 from California. We observed no structure, with an overall $F_{st}$
#'   value of -0.021. We conclude that pea aphid populations in the United
#'   States function as a single, panmictic population." 
#' * "[$\theta_w$ and $\theta_{\pi}$] for all sites across the genome were 0.0050 and 
#'   0.0045, respectively"
#' 
#' 
#' # Initial information
#' 
#' 
#' From the paper's information above, we have...
#' 
theta_w <- 0.0050
theta_pi <- 0.0045
#' 
#' 
#' I'm going to simulate a sample size of 10.
#' (The period is prepended to avoid conflicts with other object names.)
#' 
.n <- 10
#' 
#' 
#' 
#' 
#' # Calculating parameters from paper data
#' 
#' The two main pieces of information I want for the calculations are 
#' (1) the proportion of segregating sites and
#' (2) some measure of how different individuals are at segregating sites.
#' 
#' ## Segregating sites
#' 
#' [Watterson's (1975)](http://linkinghub.elsevier.com/retrieve/pii/0040580975900209)
#' estimator ($\theta_w$) is as follows:
#' 
#' $$\theta_w = \frac{ K }{ a_n }$$
#' 
#' where $K$ is the proportion of segregating sites. Variable $a_n$ is below:
#' 
#' $$a_n = \sum_{i=1}^{n-1} \frac{1}{i}$$
#' 
#' where $n$ is the number of individuals sampled.
#' Thus the the proportion of segregating sites is simply $K = \theta_w a_n$.
#' 
#' For simplicity, I'm going to first make a function to compute $a_n$ for a given 
#' value or values of $n$.
#' 
a_n <- function(n) {
    # Inner function to get a single "harmonic number"
    harm_n <- function(inner_n) {
        harm_n_vec <- 1 / 1:inner_n
        return(sum(harm_n_vec))
    }
    if (any((n %% 1) != 0)) stop("n must be entirely integers")
    sapply(n - 1, harm_n)
}
#' 
#' 
#' Now to compute the proportion of segregating sites for my chosen sample size of 
#' `r .n`, I simply multiply $\theta_w$ and $a_{10}$.
#' 
(seg_sites <- round(theta_w * a_n(10), 5))
#' 
#' 
#' 
#' 
#' 
#' ## Diversity at segregating sites
#' 
#' 
#' [Nei and Li's (1979)](http://www.pnas.org/cgi/doi/10.1073/pnas.76.10.5269)
#' measure of nucleotide diversity, $\theta_{\pi}$, is calculated using the following
#' equation:
#' 
#' $$\theta_\pi = \sum_{ij} x_i x_j \pi_{ij}$$
#' 
#' Here, $x_i$ and $x_j$ represent the frequencies of the $i$th and $j$th unique 
#' sequences respectively and 
#' $\pi_{ij}$ represents the proportion of divergent sequence between the $i$th and 
#' $j$th unique sequences.
#' 
#' If I assume that all lines will be unique sequences—a safe assumption if whole 
#' genomes are considered—then the above equation can be expressed as follows:
#' 
#' $$\theta_\pi = \sum_{ij} \frac{1}{n^2} \pi_{ij}$$
#' 
#' Then, since the number of total pairwise combinations between $n$ sequences can
#' be simplified to ${n \choose 2}$, we can use $\bar{\pi}$, the mean proportional 
#' sequence divergence between any two sequences, as such:
#' 
#' $$\theta_\pi = {n \choose 2} \frac{1}{n^2} \bar{\pi}$$
#' 
#' Solving for $\bar{\pi}$ yields the following:
#' 
#' $$\bar{\pi} = \frac{\theta_\pi n^2}{{n \choose 2}}$$
#' 
#' 
#' Since I've already calculated the proportion of segregated sites, I want the mean
#' divergence at segregated sites only.
#' (This improves computational and coding efficiency because I only have to worry about
#' segregating sites.)
#' To only consider segregated sites, I divide the whole expression by the proportion 
#' of segregated sites.
#' This leaves me with the proportional nucleotide divergence between two sequences
#' at segregating sites.
#' 
(seg_div <- round({theta_pi * .n^2 / choose(.n, 2)} / seg_sites, 4))
#' 
#' 
#' 
#' 
#' # Functions to create 
#' 
#' 
#' These functions (1) do pairwise comparisons for a vector of sequences (each sequence
#' containing one nucleotide for each sample) and (2) create a sequence from a 
#' specified number of each nucleotide.
#' 
pw_comp <- function(seq_vector) {
    if (length(seq_vector[[1]]) == 1) {
        seq_vec <- .Internal(strsplit(seq_vector, '', FALSE, FALSE, FALSE))
    }
    output <- sapply(seq_vec, 
                     function(.s) {
                         pw_mat <- t(combn(.s, 2))
                         diffs <- ifelse(pw_mat[,1] == pw_mat[,2], 0, 1)
                         return(mean(diffs))
                     })
    return(output)
}
make_seq <- function(a, c, g, t, return_vec = FALSE) {
    seq_vec <- c(rep('A', as.integer(a)), rep('C', as.integer(c)), 
                 rep('G', as.integer(g)), rep('T', as.integer(t)))
    if (return_vec) return(seq_vec)
    return(paste(seq_vec, collapse = ''))
}
#' 
#' 
#' This function creates a nucleotide frequency matrix containing, in each row, 
#' nucleotide frequencies that (1) coincides closest with
#' specified divergence at segregated sites and (2) sum to the number of samples.
#' The order of the frequencies is not considered important, since I intend
#' rows in this frequency matrix to be shuffled when used.
#' 
#' This takes ~1 min for `N=100`, ~2 sec for `N=50`, and << 1 sec for `N=10`.
#' 
nt_freq <- function(N, divergence) {
    freq_mat <- combinations(N + 1, 4, 0:N, set = FALSE, repeats.allowed = TRUE)
    freq_sums <- rowSums(freq_mat)
    freq_mat <- freq_mat[freq_sums == N,]
    pw_divs <- pw_comp(apply(freq_mat, 1, function(x) do.call(make_seq, as.list(x))))
    min_diff_inds <- which(abs(pw_divs - divergence) == min(abs(pw_divs - divergence)))
    cat(sprintf('Minimum absolute difference from divergence = %s\n', 
                format(min(abs(pw_divs - divergence)), scientific = TRUE, digits = 4)))
    return(freq_mat[min_diff_inds,])
}
#' 
#' 
#' 
#' 
#' This function returns a `seq_obj` object containing an innner character vector,
#' matrix, and numeric value.
#' The first 2 inner objects contain information on each input sequence, and positions in 
#' output objects coincide with positions in the input `DNAStringSet` object.
#' So information in row `i` of `freq_len` and position `i` in `seqs` contains info
#' for the `i`th sequence in the input `DNAStringSet` object.
#' The matrix `freq_len` contains two columns, one for the number of segregating sites
#' and one for the length of the sequence.
#' The character vector `seqs` contains the sequences themselves.
#' The numeric value (`N`) contains the total number of sequences.
#' 
#' These `seq_obj` objects are used for downstream processes, and this function is not
#' intended to be used outside another function.
#' 
#+ make_constr_objs
setClass(Class = 'seq_obj', representation(seqs = 'character', freq_len = 'matrix', 
                                           N = 'numeric'))
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
    seq_obj <- new('seq_obj', freq_len = freq_len, seqs = seqs, N = length(seqs))
}
#' 
#' The inner objects can be accessed as such, for a `seq_obj` named `Z`:
#' 
#' ```{r access_seq_obj, eval = FALSE}
#' Z@freq_len
#' Z@seqs
#' Z@N
#' ```
#' 
#' 
#' 
#' This function randomly choses locations within a sequence for a given row in a 
#' `seq_obj@freq_len` matrix.
#' 
one_sites <- function(.freq_len_row) {
    samp_n <- .freq_len_row[1]
    if (samp_n == 0) {
        return(integer(0))
    }
    seq_length <- .freq_len_row[2]
    out_locs <- .Internal(sample(seq_length, samp_n, FALSE, NULL))
    return(out_locs)
}
#' 
#' 
#' This function changes one sequence to `N` sequences, where `N` is the number of 
#' samples.
#' It also introduces variants at the supplied positions.
#' The variants will conform to the nucleotide frequency matrix that should be created
#' beforehand (see above for its description).
#' Because this function will be run many times, I wrote it in C++ (function 
#' `cpp_change` in file `variants.cpp`).
#' The R version is kept here to show the general process.
#' 
change_sites <- function(seq, positions, freq_mat, n_samps) {
    
    if (length(positions) == 0) return(rep(seq, n_samps))
    
    seq_vec <- unlist(.Internal(strsplit(seq, '', FALSE, FALSE, FALSE)))
    
    seq_mat <- matrix(rep(seq_vec, n_samps), nrow = n_samps, byrow = TRUE)
    
    ran_freq <- freq_mat[sample(nrow(freq_mat), 1), sample(4)]
    new_nucs <- sapply(positions, 
                       function(i) do.call(make_seq, as.list(c(ran_freq, TRUE))))
    
    seq_mat[,positions] <- new_nucs
    
    seq_out <- apply(seq_mat, 1, paste0, collapse = '')
    
    return(seq_out)
}


#' 
#' # Testing functions
#' 
#' I am `source`-ing `wr_size_filter.R` to use those objects to first filter the 
#' fragments by size before removing faraway sequences.
#' 
#+ trick_preamble, echo = FALSE
# This is to keep `wr_size_filter.R` from source-ing `wr_preamble.R`
.preamble_sourced <- TRUE
#' 
#+ source_size_filter
source('wr_size_filter.R')
#' 
dna_ss <- sread(readFasta(sprintf('./genome_data/frags_%s.fa.gz', 'ApeKI'))) %>% 
    size_filter
divergence = seg_div
seg_prop = seg_sites
n_samps = 10




set.seed(9)
freq_mat <- nt_freq(n_samps, divergence)

seq_obj <- constr_objs(dna_ss, seg_prop)
seq_obj@seqs %>% head
seq_obj@freq_len %>% head

process_one <- function(i) {
    freq_len_row <- seq_obj@freq_len[i,]
    seq <- seq_obj@seqs[i]
    sites <- one_sites(freq_len_row)
    new_seqs <- cpp_change(seq, sites-1, freq_mat, n_samps)
    return(new_seqs)
}
new_sites <- mclapply(1:seq_obj@N, process_one, mc.cores = 3)
new_sites <- c(new_sites, recursive = TRUE)



# sourceCpp('variants.cpp')



system.time({new_sites <- mclapply(1:seq_obj@N, process_one, mc.cores = 3) %>% 
    c(recursive = TRUE)})








# cpp_change(seq, sites-1, freq_mat, n_samps)

# cpp_par_newseqs(freq_len, seqs, freq_mat, n_samps)






# system.time({rand_locs <- get_sites(freq_len)})
# system.time({rand_locs <- cpp_get_sites(freq_len)})

pos_mat <- cbind(rand_seqs, rand_locs)
head(pos_mat)

length(rand_locs3)



change_all_sites <- function(seq_vec, in_seq_nums, positions_mat, freq_mat) {
    seq_nums <- positions_mat[,1]
    positions <- positions_mat[,2]
    n_samps <- sum(freq_mat[1,])
    if (length(in_seq_nums) != length(seq_vec)) {
        stop('in_seq_nums should be same length as seq_vec')
    }
    out_list <- 
        mclapply(1:length(seq_vec),
                 function(.i) {
                     seq <- seq_vec[.i]
                     seq_n <- in_seq_nums[.i]
                     positions_i <- positions[seq_nums == seq_n] - 1
                     new_seqs <- cpp_change(seq, positions_i, freq_mat, n_samps)
                     return(new_seqs)
                 }, mc.cores = 3)
    out_vec <- c(out_list, recursive = TRUE)
    return(out_vec)
}



sub <- round(length(char_fasta) / 10)
set.seed(8); sub_seqs <- sample(length(char_fasta), sub)
sub_fasta <- char_fasta[sub_seqs]
sub_pos <- pos_mat[pos_mat[,1] %in% sub_seqs,]

set.seed(1); system.time({r_test <- change_all_sites(sub_fasta, sub_seqs, sub_pos, freq_mat)})
#    user  system elapsed
#   2.109   0.918   0.975














# /////////////////////////////////////
# /////////////////////////////////////

# Left off: Putting the above code together cogently and efficiently

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# // #include <algorithm>

library(Rcpp)
sourceCpp(code = 
'// #include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]


')




# This is an older C++ version of the function in dv_make_strands.R
# library(Rcpp)
# cppFunction(
#     'std::vector<std::string> cut_seqs(std::vector<std::string> seq, int cut_len) {
#         int n = seq.size();
#         int seq_len;
#         std::vector<std::string> output(2 * n);
#         for(int i = 0, j = 0; i < n; i++, j+=2) {
#             seq_len = seq[i].length();
#             if(seq_len < cut_len){
#                 output[j] = seq[i];
#             } else {
#                 output[j] = seq[i].substr(0, cut_len);
#                 output[j+1] = seq[i].substr(seq_len - cut_len, cut_len);
#             }
#         }
#         return output;
#     }'
# )



# compare_pw <- function(cpp_test, r_test) {
#     pw_df <- lapply(
#         seq(1, length(cpp_test)-9, 10), 
#         function(i) {
#             cpp_i <- cpp_test[i:(i+9)] %>% 
#                 lapply(function(.s) unlist(.Internal(strsplit(.s, '', F, F, F)))) %>% 
#                 do.call(what = cbind, args = .) %>% 
#                 apply(1, cpp_merge_str) %>% 
#                 pw_comp %>% mean
#             r_i <- r_test[i:(i+9)] %>% 
#                 lapply(function(.s) unlist(.Internal(strsplit(.s, '', F, F, F)))) %>% 
#                 do.call(what = cbind, args = .) %>% 
#                 apply(1, cpp_merge_str) %>% 
#                 pw_comp %>% mean
#             data.frame(cpp = cpp_i, r = r_i)
#         }) %>% 
#         bind_rows %>% 
#         as.tbl %>% 
#         mutate(diff = cpp - r)
#     return(pw_df)
# }


#' <!--- References for reading and writing fastas
# dig_frags <- lapply(setNames(.chosen_enz, .chosen_enz), function(enz) {
#     fasta <- readFasta(sprintf('./genome_data/frags_%s.fa.gz', enz))
#     sread(fasta)
# })
# writeFasta(write_list[[enz]], file = sprintf('./genome_data/frags_%s.fa.gz', enz), 
#            mode = 'w', compress = 'gzip')
# cat(sprintf('%s file finished', enz), '\n')
#' -->





#' (See the `README.md` file for why I'm including `./genome_data/` in file paths.)










