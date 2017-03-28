# This preamble loads packages and runs code necessary for >1 "working R" files

suppressPackageStartupMessages({
    
})
library(SimRAD)
library(dplyr)
library(purrr)
library(ggplot2)

# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))

