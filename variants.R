#' ---
#' title: 'Create structural variants from reference genome'
#' author: 'Lucas Nell'
#' date: '24 March 2017'
#' output: github_document
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
#' To install `RSVSim`, run the following code:
#' 
#' ```{r, eval = FALSE}
#' source('https://bioconductor.org/biocLite.R')
#' biocLite('RSVSim')
#' ```
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
        library(RSVSim)
    })
#' 
#' 
#' 



#' (See the `README.md` file for why I'm including `./genome_data/` in file paths.)











