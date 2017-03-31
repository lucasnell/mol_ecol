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
})
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
#' This function creates a nucleotide frequency matrix that coincides closest with
#' specified divergence at segregated sites and sums to the number of samples.
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
#' This function changes one sequence to `N` sequences, where `N` is the number of 
#' samples.
#' It also introduces variants at the supplied positions.
#' The variants will conform to the nucleotide frequency matrix that should be created
#' beforehand (see above for its description).
#' This took 0.012 sec to change 50 positions in a 447-length sequence.
#' 
#' 
change_sites <- function(seq, positions, freq_mat) {
    
    N <- sum(freq_mat[1,])
    
    if (length(positions) == 0) return(rep(seq, N))
    
    seq_vec <- unlist(.Internal(strsplit(seq, '', FALSE, FALSE, FALSE)))
    
    seq_mat <- matrix(rep(seq_vec, N), nrow = N, byrow = TRUE)
    
    ran_freq <- freq_mat[sample(nrow(freq_mat), 1), sample(4)]
    new_nucs <- sapply(positions, 
                       function(i) do.call(make_seq, as.list(c(ran_freq, TRUE))))
    
    seq_mat[,positions] <- new_nucs
    
    seq_out <- apply(seq_mat, 1, paste0, collapse = '')
    
    return(seq_out)
}
system.time(replicate(1000, change_sites(char_fasta[1], seq.int(1, 130, 2), freq_mat)))
system.time(replicate(1000, cpp_change_sites(char_fasta[1], seq.int(1, 130, 2), freq_mat)))
#' 
#' 
#' 
#' 
#' 
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
test_fasta <- sread(readFasta(sprintf('./genome_data/frags_%s.fa.gz', 'ApeKI'))) %>% 
    size_filter
char_fasta <- as.character(test_fasta)
freq_mat <- nt_freq(10, seg_div)

# get_seg_sites <- function(seq_lens, seg_sites)
# }

library(Rcpp)
library(RcppArmadillo)
sourceCpp('variants.cpp')


get_sites <- function(.seq_num, .samp_n, .seq_lens) {
    ran_locs <- .Internal(sample(.seq_lens[.seq_num], .samp_n, FALSE, NULL))
    return(ran_locs)
}

seq_lens <- nchar(char_fasta)
total_seg <- round(sum(seq_lens) * seg_sites, 0)
set.seed(8)
rand_seqs <- sort(sample(length(seq_lens), total_seg, replace = TRUE, prob = seq_lens))
seq_freq <- mutate(as_data_frame(table(rand_seqs)), seq = as.integer(rand_seqs)) %>% 
    select(seq, n) %>% 
    as.matrix
head(seq_freq)
system.time({rand_locs <- c(mapply(get_sites, seq_freq[,1], seq_freq[,2],
                                   MoreArgs = list(.seq_lens = seq_lens)), 
                            recursive = TRUE)})
system.time({rand_locs <- cpp_get_sites(seq_lens, seq_freq)})

pos_mat <- cbind(rand_seqs, rand_locs)
head(pos_mat)



sub <- round(length(char_fasta) / 10000)
sub_fasta <- char_fasta[1:sub]
sub_pos <- pos_mat[pos_mat[,1] <= sub,]



change_all_sites <- function(seq_vec, positions_mat, freq_mat) {
    seq_nums <- positions_mat[,1]
    positions <- positions_mat[,2]
    out_list <- lapply(1:length(seq_vec), 
                       function(.i) {
                           seq <- seq_vec[.i]
                           positions_i <- positions[seq_nums == .i]
                           new_seq <- change_sites(seq, positions_i, freq_mat)
                           return(new_seq)
                       })
    out_vec <- c(out_list, recursive = TRUE)
    return(out_vec)
}

system.time({set.seed(1); cpp_test <- cpp_change_sites(sub_fasta, sub_pos, freq_mat)})
system.time({set.seed(1); r_test <- change_all_sites(sub_fasta, sub_pos, freq_mat)})
identical(cpp_test, r_test)













# /////////////////////////////////////
# /////////////////////////////////////

# Left off: Putting the above code together cogently and efficiently

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\




sourceCpp(code = 
'#include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]]
IntegerVector match__( IntegerVector x, int i){
    IntegerVector y;
    y = x[x == i];
    return y.size();
}
')


match__(1:5, 6)







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










