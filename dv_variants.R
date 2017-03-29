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
(seg_sites <- theta_w * a_n(10))
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
(seg_div <- (theta_pi * .n^2 / choose(.n, 2)) / seg_sites)
#' 
#' 
#' 
#' 
#' # Functions to create 
#' 

fasta <- sread(readFasta(sprintf('./genome_data/frags_%s.fa.gz', 'BspEI')))
# fasta[[1]] <- DNAString('AA')
# DNAStringSet(c('AA', 'GG'))

# get_seg_sites <- function(seq_lens, seg_sites)
# }

seq_num = 1
.one_seq <- function(seq_num) {
    system.time({seq_length <- seq_lens[seq_num]
    samp_n <- length(rand_seqs[rand_seqs == seq_num])
    ran_locs <- sample(1:seq_length, samp_n, replace = FALSE)})
    return(ran_locs)
}

(0.057 * 40e3) / 60


seq_lens <- width(fasta)
total_seg <- round(sum(seq_lens) * seg_sites, 0)
rand_seqs <- sort(sample(1:length(seq_lens), total_seg, replace = TRUE, prob = seq_lens))
# rand_locs <- sapply(1:length(seq_lens), .one_seq)
rand_locs <- lapply(1:40, .one_seq)


# /////////////////////////////////////
# /////////////////////////////////////

# Left off: Above took a long time, so I moved onto filtering out areas not adjacent to 
#           cut sites since this will reduce the fasta file sizes.

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



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










