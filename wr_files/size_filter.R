
# These are objects necessary for filtering sequences based on sequence size.

# If wr_preamble.R file hasn't been sourced, it needs to be.
# wr_preamble.R creates the .preamble_sourced object when it's run
if(!'.preamble_sourced' %in% ls(all.names = TRUE)) source('./wr_files/preamble.R')

size_filter <- function(dna_ss) {
    .prob_coef <- 594.6
    .prob_dens <- function(.frag_sizes) dnbinom(.frag_sizes, size = 2.013, mu = 144.8)
    .fast_bern <- function(.p) {
        .U <- runif(length(.p),0,1)
        .outcomes <- .U < .p
        return(.outcomes)
    }
    frag_sizes <- width(dna_ss)
    frag_probs <- .prob_dens(frag_sizes) * .prob_coef
    frag_keep <- .fast_bern(frag_probs)
    return(dna_ss[frag_keep])
}


