# This preamble loads packages and runs code necessary for >1 "working R" files

suppressPackageStartupMessages({
    library(ShortRead)
    library(magrittr)
    library(dplyr)
    library(purrr)
    library(ggplot2)
    library(Rcpp)
})

# This is for other files to check if this file has been sourced
.preamble_sourced <- TRUE

# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))

