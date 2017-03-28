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
#' This script provides background for the function that removes regions from digested
#' fragments that will not be sequenced, those that are far from restriction enzyme cut
#' sites 
#' 
#' 
#' 
#' Note: You only have to account for barcode and enzyme overhang in subtracting from 
#' read length (e.g., for 100bp reads in ApeKI using 4-bp barcodes, it's 100 - 4 - 3).
#' 
#' 
#' 
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages({
    library(purrr)
    library(dplyr)
    library(ShortRead)
})
#' 
#' 
#' 
#' 
#' 
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


