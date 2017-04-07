
# These are objects necessary for filtering sequences based on sequence size.

# If wr_preamble.R file hasn't been sourced, it needs to be.
# wr_preamble.R creates the .wr_env object when it's run
if(!'.wr_env' %in% ls(all.names = TRUE)) source('./wr_files/preamble.R')

size_filter <- function(dna_in) {
    .prob_coef <- 594.6
    .prob_dens <- function(.frag_sizes) dnbinom(.frag_sizes, size = 2.013, mu = 144.8)
    .fast_bern <- function(.p) {
        .U <- runif(length(.p),0,1)
        .outcomes <- .U < .p
        return(.outcomes)
    }
    frag_sizes <- nchar(dna_in)
    frag_probs <- .prob_dens(frag_sizes) * .prob_coef
    frag_keep <- .fast_bern(frag_probs)
    return(dna_in[frag_keep])
}


