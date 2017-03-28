# This preamble loads packages and runs code necessary for >1 "working R" files

suppressPackageStartupMessages({
    library(SimRAD)
    library(dplyr)
    library(purrr)
    library(ggplot2)
})

# This is for other files to check if this file has been sourced
.preamble_sourced <- TRUE

# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))

# Rounds proportion data to 4 digits while having the proportions still sum to 1.
rounded_p <- function(vec, round_digits = 4) {
    new_vec <- round(vec / sum(vec), round_digits)
    # If the rounded numbers don't sum to 1, then add the difference to the 
    # rounded number that differs most from the original, in the same direction as the
    # difference between sum(new_vec) and 1
    if (sum(new_vec) != 1) {
        newish_vec <- vec / sum(vec)
        diffs <- (newish_vec - new_vec) * 10^round_digits
        ind <- which(diffs == ifelse(sum(new_vec) < 1, min(diffs), max(diffs)))
        new_vec[ind] <- new_vec[ind] + (1 - sum(new_vec))
    }
    return(new_vec)
}
