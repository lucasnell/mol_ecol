#' ---
#' title: 'Filter digested fragments'
#' author: 'Lucas Nell'
#' date: '19 March 2017'
#' output: github_document
#' ---
#' 
#' *Updated `r format(Sys.Date(), '%d %B %Y')`*
#' 


#' 
#' > From six sequencing lanes, we identified 809,651 sequence tags (at least five times)
#' > from one or both flanks of 654,998 of the 2.1 million ApeKI cut sites lying within
#' > the single copy genomic fraction.
#' 
#' (Elshire et al. 2011, p 5)
#' 
#' From above, I've created below an object storing the proportion of cut sites that 
#' I'll assume get sequenced:
#' 

seq_p <- 654998 / 2.1e6




#' 
#' Used WebPlotDigitizer (arohatgi.info/WebPlotDigitizer/) to extract data from Figure 3 
#' in Elshire et al. (2011).
#' 

rounded_p <- function(vec, round_digits = 4) {
    new_vec <- round(vec / sum(vec), round_digits)
    if (sum(new_vec) != 1) {
        newish_vec <- vec / sum(vec)
        diffs <- (newish_vec - new_vec) * 10^round_digits
        ind <- which(diffs == ifelse(sum(new_vec) < 1, min(diffs), max(diffs)))
        new_vec[ind] <- new_vec[ind] + (1 - sum(new_vec))
    }
    return(new_vec)
}

frag_sizes <- read_csv('frag_sizes.csv', col_types = 'cd') %>% 
    mutate(size = as.integer(rep(seq(50, 1000, 50), each = 2)),
           type = rep(c('predicted', 'actual'), 1000 / 50)) %>% 
    select(size, type, prop) %>% 
    group_by(type) %>%
    mutate(prop = rounded_p(prop)) %>%
    ungroup


frag_sizes %>% 
    ggplot(aes(size, prop)) + 
    geom_bar(aes(fill = type), stat = 'identity', position = 'dodge') +
    theme_classic()

frag_sizes %>% 
    filter(type == 'actual') %>% 
    ggplot(aes(size, prop)) + 
    geom_bar(stat = 'identity', position = 'dodge', fill = 'dodgerblue') +
    theme_classic()


# To fit distribution to censored data:
# ?fitdistcens
# https://cran.r-project.org/web/packages/fitdistrplus/fitdistrplus.pdf



frag_sizes %>% 
    filter(type == 'actual') %>% 
    mutate(size = ifelse(size == 1000, NA, size))






#' 
#' 
#' # References
#' 
#' Elshire, R. J., J. C. Glaubitz, Q. Sun, J. A. Poland, K. Kawamoto, E. S. Buckler, 
#' and S. E. Mitchell. 2011. A robust, simple genotyping-by-sequencing (GBS) approach 
#' for high diversity species. *PLOS ONE* __6__:e19379.
#' 
#' 
