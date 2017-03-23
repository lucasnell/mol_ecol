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
#' This script accounts for the fact that in GBS, many potential cut sites are not 
#' sequenced, and that the disparity is largely driven by fragment size. 
#' It filters the digested fragments from the previous script 
#' (`digest_genome.R`).
#' 
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages(
    suppressWarnings({
        library(fitdistrplus)
        library(tidyr)
        library(purrr)
        library(readr)
        library(ggplot2)
        library(dplyr)
        library(magrittr)
        library(ShortRead)
    }))
#' 
#' 
#' 
#' 
#' # Characterizing sequenced fragments
#' 
#' 
#' ## Proportion sequenced
#' 
#' > From six sequencing lanes, we identified 809,651 sequence tags (at least five times)
#' > from one or both flanks of 654,998 of the 2.1 million *ApeKI* cut sites lying within
#' > the single copy genomic fraction.
#' 
#' ([Elshire et al. 2011](http://dx.plos.org/10.1371/journal.pone.0019379), p 5)
#' 
#' From above, I've created below an object storing the proportion of cut sites that 
#' I'll assume get sequenced when using *ApeKI* in maize 
#' (as was done in the above study):
#' 
#' 
.seq_p <- round(654998 / 2.1e6, 4)
#' 
#' 
#' 
#+ rounded_p_fun, echo = FALSE
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
#' 
#' 
#' ## Data on distribution of fragments
#' 
#' I used WebPlotDigitizer (arohatgi.info/WebPlotDigitizer/) to extract data 
#' from Figure 3 in Elshire et al. (2011) to a csv file (`frag_sizes.csv`).
#' 
#' The below code cleans up the csv file. The `rounded_p` function rounds proportion data
#' to 4 digits while having the summed proportions still add up to 1.
#' 
#+ frag_sizes_df
frag_sizes <- read_csv('frag_sizes.csv', col_types = 'cd') %>% 
    mutate(size = as.integer(rep(seq(50, 1000, 50), each = 2)),
           type = rep(c('predicted', 'sequenced'), 1000 / 50)) %>% 
    select(size, type, prop) %>% 
    group_by(type) %>%
    mutate(prop = rounded_p(prop)) %>%
    ungroup
#' 
#' 
#' This figure showed the difference between the fragment-size distribution 
#' (1) predicted from cut site locations in the maize genome and 
#' (2) from locations where reads were actually found when sequencing was performed.
#' The figure is reproduced below.
#' 
#+ frag_sizes_plot, echo = FALSE
frag_sizes %>% 
    mutate(type = factor(type, levels = c('predicted', 'sequenced'))) %>% 
    ggplot(aes(size - 25, prop)) + 
    geom_bar(aes(fill = type), stat = 'identity', width = 45,
             position = position_dodge(width = 30)) +
    theme_classic() +
    scale_fill_manual(NULL, values = c('dodgerblue', 'red')) +
    theme(legend.position = c(0.75, 0.75)) +
    scale_x_continuous('Fragment size (bp)', breaks = seq(0, 1000, 50),
                       labels = c(seq(0, 950, 50), 'Inf')) +
    ylab('Proportion of total fragments') +
    labs(caption = expression(italic(Note) * ': x-axis represents bin bounds'))
#' 
#' 
#' 
#' 
#' ## Fitting a distribution for sequenced fragments
#' 
#' I next want to fit a negative binomial distribution to the size-distribution of
#' sequenced fragments.
#' Because these data are binned (i.e., censored), I have to use the `fitdistcens`
#' function from the `fitdistrplus` package.
#' This function requires as input a data frame with columns named "left" and "right".
#' 
#' The code below creates the left and right values appropriately for each size bin
#' (see `?fitdistcens`), then creates a data frame with $p \times N$ rows of the
#' left and right values, where $N$ is the total number of fragments and $p$ is the
#' proportion of total fragments represented by that size bin (i.e., the height
#' of its bar in the plot above).
#' I'm assuming here that the total number of fragments ($N$) is the 809,651 from the 
#' Elshire et al. (2011) quote at the top of this document.
#' 
cens_df <- frag_sizes %>% 
    filter(type == 'sequenced') %>% 
    split(.$size) %>% 
    map_df(
        function(.x, total_frags = 809651) {
            .left <- .x$size - 49
            .right <- ifelse(.x$size == 1000, NA, .x$size)
            N <- round(.x$prop * total_frags, 0)
            data_frame(left = rep(.left, N),  right = rep(.right, N))
        }
    )
#' 
#' 
#' Now I do the actual fitting using `fitdistcens`. (This took ~1 minute)
#' 
#' *Note*: Even though `fitdistcens` can't use a "tibble" data frame, I'm keeping
#' `cens_df` in that format to keep me from accidentally printing >800,000 rows.
#' 
#' 
seq_fit <- fitdistcens(as.data.frame(cens_df), 'nbinom')
summary(seq_fit)
#' 
#' 
#' 
#' I can now create a function to calculate the probability density for a given 
#' fragment size.
#' 
.prob_dens <- function(frag_sizes) {
    dnbinom(frag_sizes, size = 2.013, mu = 144.8)
}
#' 
#' 
#' Here is the probability density function (black curve, right y-axis) along with
#' the binned distribution of fragment sizes for sequenced reads 
#' (blue bars, left y-axis):
#' 
#+ frag_size_plot, echo = FALSE
frag_sizes %>% 
    filter(type == 'sequenced') %>% 
    ggplot(aes(y = prop)) + 
    geom_bar(aes(size - 25), stat = 'identity', position = 'dodge', 
             fill = 'dodgerblue') +
    theme_classic() + 
    geom_line(data = data_frame(size = 1:1000, prop = .prob_dens(size)) %>% 
                  mutate(prop = max(frag_sizes$prop) * prop / max(prop)),
              aes(size), size = 0.75) +
    scale_x_continuous('Fragment size (bp)', breaks = seq(0, 1000, 50),
                       labels = c(seq(0, 950, 50), 'Inf')) +
    scale_y_continuous(
        'Proportion of total fragments', 
        sec.axis = sec_axis(~ . * max(.prob_dens(1:1000)) / max(frag_sizes$prop), 
                            name = 'Density')) +
    labs(caption = expression(italic(Note) * ': x-axis represents bin bounds for bars'))
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Reading in digestion fragments
#' 
#' 
#' I now read the fasta files made in `digest_genome.R` that represent in silico 
#' digestions of the aphid genome using three different restriction enzymes.
#' 
#' 
chosen_enz <- c('ApeKI', 'BstBI', 'NruI-HF')
dig_frags <- lapply(setNames(chosen_enz, chosen_enz), function(enz) {
    fasta <- readFasta(sprintf('frags_%s.fa.gz', enz))
    sread(fasta)
})
#' 
#' 
#' This is a data frame containing, for each restriction enzyme, all fragment lengths and
#' their probability densities from the sequenced-fragment probability density function.
#' 
dig_frag_df <- lapply(chosen_enz, 
                      function(re){
                          .frag_len <- width(dig_frags[[re]])
                          .prob <- .prob_dens(.frag_len)
                          .frag_index <- 1:length(.frag_len)
                          .df <- data_frame(enzyme = re, 
                                            frag_index = .frag_index,
                                            frag_len = .frag_len, 
                                            prob = .prob)
                          return(.df)
                      }) %>% 
    bind_rows
#' 
#' 
#' 
#' 
#' 
#' # Calibrating probability distribution for *ApeKI*
#' 
#' 
#' From here, I could filter *ApeKI*-digested fragments using the following code:
#' 
#' ```{r sample_code, eval = FALSE}
#' dig_frag_df %>% 
#'     filter(enzyme == 'ApeKI') %>% 
#'     sample_frac(.seq_p, weight = prob)
#' ```
#' 
#' ... which would give fragments with a higher probability density a greater chance
#' of being selected and would give me the proportion I expect based on the maize-genome
#' data.
#' 
#' If I did that for the other restriction enzymes, I would have to assume that the same
#' proportion of fragments would be sequenced with them as with *ApeKI*.
#' These restriction enzymes created longer average fragments than *ApeKI*, and we're
#' assuming that fragment size is the primary determinant of whether a fragment is
#' sequenced. Thus I decided not to use the sampling approach above.
#' 
#' What I decided upon was to multiply the probability densities by a coefficient while
#' limiting the new resulting probability to a maximum of 1:
#' 
#' $$p_i^{\prime} = \min(\{ \beta p_i, 1 \})$$
#' 
#' where 
#' $\beta$ is the selected coefficient and
#' $p_i$ is the probability density of fragment $i$.
#' I treated whether fragment $i$ is sequenced as a Bernoulli trial with probability
#' $p_i^{\prime}$.
#' 
#' The expected proportion selected ($\mathbb{E}(P)$) was calculated as such:
#' 
#' $$\mathbb{E}(P) = \frac{\sum_{i=1}^{n} p_i^{\prime} }{n}$$
#' 
#' where $n$ is the total number of fragments.
#' 
#' I tested multiple coefficients to find the one that minimized the absolute difference
#' between the proportion of individuals we would predict would be selected from our 
#' data (i.e., $\mathbb{E}(P)$) and the proportion sequenced from the maize data ($P_m$).
#' 
#' 
#' 
#' 
#+ get_multiplier, echo = FALSE
# (Since these objects are temporary, I'm making them start with 'test_' and end in '_' so
# they're easy to search for and remove.)
test_seq_ <- seq(594.5, 595, 0.1)
test_prob_ <- filter(dig_frag_df, enzyme == 'ApeKI')$prob
test_m_ <- sapply(test_seq_,
                  function(m) abs(mean(ifelse(test_prob_ * m > 1, 1, test_prob_ * m)) - 
                                      .seq_p))
par(mar = c(5, 4.5, 1, 1))
plot(test_m_ ~ test_seq_, type = 'l', ylab = expression('|' ~ E(P) - P[m] ~ '|'), 
     xlab = 'Multiplier')
.prob_coef <- test_seq_[test_m_ == min(test_m_)]
rm(list = ls()[grepl('^test_.*_$', ls())])
#' 
#' 
#' The best coefficient (`r .prob_coef`) was assigned to object `.prob_coef`.
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Testing fragment filtering by size
#' 
#' 
#' In this section I test how the calibrated probability densities impact the 
#' proportion of fragments sequenced for all enzymes, using simulated Bernoulli trials.
#' 
#' The `.fast_bern` function quickly performs many Bernoulli trials, each with a unique 
#' probability (link to where I found it
#' [here](http://r-bloggers.com/variable-probability-bernoulli-outcomes-fast-and-slow/)).
#' It is ~75 faster than using `sapply` on a `rbinom(1,1,p)` call, and it returns a 
#' logical instead of a binary integer vector, which makes filtering simpler.
#' 
#' The `test_filter` function filters the `dig_frag_df` data frame by values in the
#' `prob` column multiplied by `.prob_coef`, then it gets the number of rows present in 
#' the filtered data frame.
#' It returns the proportion of fragments retained for sequencing.
#' 
.fast_bern <- function(p) {
    U <- runif(length(p),0,1)
    outcomes <- U < p
    return(outcomes)
}
test_filter <- function(enz) {
    filt_df <- filter(dig_frag_df, enzyme == enz)
    N <- nrow(filter(filt_df, .fast_bern(prob * .prob_coef)))
    return(N / nrow(filt_df))
}
#' 
#' 
#' Now I run the tests, and find that, as I intended, the proportion of sequenced reads
#' decreases with restriction enzymes having fewer cut sites.
#' 
set.seed(329)
sapply(chosen_enz, 
                    function(e) sapply(1:100, function(i) test_filter(e))) %>% 
    as_data_frame %>% 
    gather(enzyme, prop) %>% 
    group_by(enzyme) %>% 
    summarize(`mean` = mean(prop), `sd` = sd(prop)) %>% 
    as.data.frame
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Creating functions to simulate size filtering
#' 
#' 
#' 
#' 
#' 
#' 



# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================

# LEFT OFF: Organize and apply filter using .prob_coef to other digestions

# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================


# # Visualizing output:
# data_frame(size = .test$frag_len[filt_frags(.test$prob)]) %>% 
#     ggplot() +
#     theme_classic() +
#     geom_histogram(aes(size, ..density..), binwidth = 50, fill = 'dodgerblue', 
#                    closed = 'left', center = 25) +
#     geom_line(data = data_frame(size = 1:750, prob = .prob_dens(size)) %>% 
#                   mutate(prob = 0.0055 * prob / max(prob)),
#               aes(size, prob), size = 0.75)





