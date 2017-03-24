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
    })
#' 
#' 
#' 

harm_n <- function(n) {
    harm_n_vec <- 1 / 1:n
    return(sum(harm_n_vec))
}


#' 
#' * "We sequenced 21 genetically distinct lines of pea aphids..."
#' * "...we calculated F<sub>st</sub> levels across the genome, comparing 11
#'   pea aphid lines from the Northeast US (New York and Massachusetts)
#'   and 10 from California. We observed no structure, with an overall F<sub>st</sub>
#'   value of -0.021. We conclude that pea aphid populations in the United
#'   States function as a single, panmictic population."
#' * "$\theta^w$ and $\pi$ for all sites across the genome were 0.0050 and 0.0045, 
#'   respectively"
#' * "[The aphid genome they used was] ... 541,675,471 bp total length."
#' 
#' ([Bickel et al. 2013](http://www.g3journal.org/content/3/6/993), p 996)
#' 
#' 


theta_w <- 0.0050
theta_pi <- 0.0045

# Below is proportion of segregating sites in their 21 samples
# (Probably want to figure out how that'll relate to proportion with 10 samples)
theta_w * harm_n(21 - 1)

57 / harm_n(95)

# 57 segregating sites in 96 samples

#' Left off --> Back calculate something useful from the theta_pi value
#' http://people.ibest.uidaho.edu/~bree/courses/14_EBR_ME_diversity.pdf
#' https://molpopgen.github.io/libsequence/doc/html/classSequence_1_1PolySNP.html#a9c26d27c7b0aaf1faf859272871b3cf5
#' http://www.g3journal.org/content/ggg/3/6/993.full.pdf
#' 






#' (See the `README.md` file for why I'm including `./genome_data/` in file paths.)










