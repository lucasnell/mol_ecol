#' ---
#' title: 'Make strands for sequence simulation'
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
#' This script tests the function that creates strands to be used for sequencing 
#' simulation.
#' It is to be used on digested fragments.
#' It removes regions from fragments that will not be sequenced (those that are far from
#' restriction enzyme cut sites).
#' For each cut-site-adjacent area, it returns a forward and reverse strand, the latter 
#' of which is returned as its reverse complement.
#' The reverse complement is returned because I will be simulating sequencing from
#' fragments unidirectionally in the forward direction.
#' I chose to do unidirectional sequencing because this more similarly simulates how
#' both forward and reverse fragment ends can be assigned sequencing primers.
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
#' The `make_strands` removes portions of fragments that are far from restriction enzyme
#' cut sites.
#' A `DNAStringSet` object of sequence fragments is the required argument.
#' Read length (defaults to 100bp) and barcode length (defaults to 4) are optional
#' arguments.
#' It returns a `DNAStringSet` object with the new sequences as described above.
#' The returned set should have twice as many sequences as the input set of sequences
#' because `make_strands` outputs a sequence for each sequencing strand.
#' 
#+ make_function
make_strands <- function(dna_ss, read_len = 100, bc_len = 4) {
    # Length to cut fragment to is simply the read length minus the barcode length and
    # the restriction enzyme's overhang length
    cut_len <- as.integer(read_len - bc_len)
    # Extracting character vector from input DNAStringSet
    character_seqs <- as.character(dna_ss)
    # Inner function that cuts (i.e., removes) faraways from one sequence
    .one_cut <- function(.s) {
        .seq_len <- as.integer(nchar(.s))
        if (.seq_len < cut_len) {
            f_strand <- r_strand <- .s
        } else {
            f_strand <- .Internal(substr(.s, 1L, cut_len))
            r_strand <- .Internal(substr(.s, .seq_len - cut_len + 1L, .seq_len))
        }
        return(matrix(c(f_strand, r_strand), nrow = 1))
    }
    # applying that inner function to each string in the character vector
    mat_list <- lapply(character_seqs, .one_cut)
    cut_seq_char <- do.call(rbind, mat_list)
    # Creating DNAStringSet objects for forward and reverse strands.
    # Because I am doing unidirectional (forward-only) sequencing, I need to make the 
    # reverse strands reverse complements of the forward one.
    f_strands <- DNAStringSet(cut_seq_char[,1])
    r_strands <- reverseComplement(DNAStringSet(cut_seq_char[,2]))
    # Now I combine them
    cut_dna_ss <- append(f_strands, r_strands)
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
#' Below is a run of `make_strands` including the run time and output.
#' 
#+ run_function
system.time({ms_test <- make_strands(filt_fasta)})
ms_test
#' 
#' 
#' 
#' # Session info
#' 
#+ make_session_info, echo = FALSE
devtools::session_info()
