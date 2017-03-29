#' ---
#' title: 'Remove areas faraway from cut sites'
#' author: 'Lucas Nell'
#' date: '28 March 2017'
#' output: 
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#'     fig_width: 8
#'     fig_height: 6
#' ---
#' 
#' *Updated `r format(Sys.Date(), '%d %B %Y')`*
#' 
#' 
#' 
#' 
#' This script tests the function that removes regions from digested
#' fragments that will not be sequenced (those that are far from restriction enzyme cut
#' sites).
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
suppressPackageStartupMessages(library(ShortRead))
#' 
#' 
#' 
#' 
#' # References
#' 
#' These are the papers I read to develop the below function:
#' 
#' Davey, J. W., P. A. Hohenlohe, P. D. Etter, J. Q. Boone, J. M. Catchen, and M. L. 
#' Blaxter. 2011. Genome-wide genetic marker discovery and genotyping using 
#' next-generation sequencing. *Nature Reviews Genetics* __12__:499–510.
#' Available from [this link](http://dx.doi.org/10.1038/nrg3012).
#' 
#' Elshire, R. J., J. C. Glaubitz, Q. Sun, J. A. Poland, K. Kawamoto, E. S. Buckler, 
#' and S. E. Mitchell. 2011. A robust, simple genotyping-by-sequencing (GBS) approach 
#' for high diversity species. *PLOS ONE* __6__:e19379.
#' Available from 
#' [this link](http://dx.plos.org/10.1371/journal.pone.0019379).
#' 
#' 
#' 
#' # Function
#' 
#' The `rm_faraway` removes portions of fragments that are far from restriction enzyme
#' cut sites.
#' DNA sequences (as a `DNAStringSet` object) and an enzyme name are the required 
#' arguments.
#' Read length (defaults to 100bp) and barcode length (defaults to 4) are optional
#' arguments.
#' It returns a `DNAStringSet` object with the faraway fragments removed.
#' The returned set should have twice as many sequences as the input set of sequences
#' because `rm_faraway` outputs a sequence for each sequencing strand.
#' 
#+ make_function
rm_faraway <- function(dna_ss, enzyme, read_len = 100, bc_len = 4) {
    overhang <- c(3, 2, 4)[c('ApeKI', 'BstBI', 'BspEI') == enzyme]
    # Length to cut fragment to is simply the read length minus the barcode length and
    # the restriction enzyme's overhang length
    cut_len <- as.integer(read_len - bc_len - overhang)
    # Extracting character vector from input DNAStringSet
    character_seqs <- as.character(dna_ss)
    # Inner function that cuts one sequence
    .one_cut <- function(.s) {
        .seq_len <- as.integer(nchar(.s))
        if (.seq_len < cut_len) {
            return(rep(.s,2))
        } else {
            return(c(.Internal(substr(.s, 1L, cut_len)), 
                     .Internal(substr(.s, .seq_len - cut_len + 1L, .seq_len))))
        }
    }
    # applying that inner function to each string in the character vector
    cut_seq_char <- unlist(lapply(character_seqs, .one_cut))
    # Now converting back to DNAStringSet
    cut_dna_ss <- DNAStringSet(cut_seq_char)
    return(cut_dna_ss)
}
#' 
#' 
#' 
#' # Testing function
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
#' Now I read the fasta file and do the filtering by size.
#' 
#+ make_filt_fasta
test_fasta <- sread(readFasta('./genome_data/frags_BstBI.fa.gz'))
filt_fasta <- size_filter(test_fasta)
#' 
#' 
#' Below is a run of `rm_faraway` including the run time and output.
#' 
#+ run_function
system.time({rm_test <- rm_faraway(filt_fasta, 'BstBI')})
rm_test
#' 
#' 
#' 
#' # Session info
#' 
#+ make_session_info, echo = FALSE
devtools::session_info()
