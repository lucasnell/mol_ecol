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
#' 
#' 
#' 
#' 
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages(
    suppressWarnings({
        library(fitdistrplus)
        library(purrr)
        library(readr)
        library(ggplot2)
        library(dplyr)
        library(magrittr)
        library(ShortRead)
    }))

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
    theme_classic() +
    scale_fill_manual(values = c('dodgerblue', 'red'))

# Data frame incorporating censoring
cens_df <- frag_sizes %>% 
    split(interaction(.$size, .$type)) %>% 
    map_df(
        function(.x) {
            .left <- .x$size - 49
            .right <- ifelse(.x$size == 1000, NA, .x$size)
            .type = .x$type
            N <- round(.x$prop * 1e4, 0)
            data_frame(type = rep(.type, N), 
                       left = rep(.left, N),  right = rep(.right, N))
        }
    )


actual_fit <- cens_df %>% 
    filter(type == 'actual') %>% 
    select(-type) %>% 
    as.data.frame %>% 
    # fitdistcens('lnorm')
    fitdistcens('nbinom')
actual_fit %>% summary





predicted_fit <- cens_df %>% 
    filter(type == 'predicted') %>% 
    select(-type) %>% 
    as.data.frame %>% 
    # fitdistcens('lnorm')
    fitdistcens('nbinom')
predicted_fit %>% summary



act_prop <- function(sizes) {
    dnbinom(sizes, size = coef(actual_fit)[['size']], 
            mu = coef(actual_fit)[['mu']])
}

pred_prop <- function(sizes) {
    dnbinom(sizes, size = coef(predicted_fit)[['size']], 
            mu = coef(predicted_fit)[['mu']])
}



bin_prop <- function(sizes, bin_size = 50, last_max = 951, distr = 'actual') {
    
    distr_coefs <- coef(get(paste0(distr, '_fit')))
    
    size_lims <- cbind(mins = seq(1, last_max, bin_size), 
                       maxes = c(seq(bin_size, last_max - 1, bin_size), Inf))
    
    size_inds <- sapply(
        sizes,
        function(s) {
            which(size_lims[,1] <= s & size_lims[,2] >= s)
        }
    )
    
    size_lims <- tryCatch(
        size_lims[size_inds,], 
        error = function(e) {
            if (grepl('invalid subscript', e)) {
                stop(paste('You probably have zeros in your data. Since having 0-length', 
                     'fragments makes no sense, I have not allowed for that here.'))
            } else stop(e)
        }
    )
    
    bin_prop_vec <- apply(size_lims, 1, 
                       function(r) {
                           .under <- pnbinom(r[1] - 1, size = distr_coefs[['size']], 
                                             mu = distr_coefs[['mu']], lower.tail = TRUE)
                           .over <- pnbinom(r[2], size = distr_coefs[['size']], 
                                            mu = distr_coefs[['mu']], lower.tail = FALSE)
                           return(1 - .over - .under)
                       })
    
    return(bin_prop_vec)
}




frag_sizes %>% 
    filter(type == 'actual') %>% 
    ggplot(aes(y = prop)) + 
    geom_bar(aes(size - 25), stat = 'identity', position = 'dodge', 
             fill = 'dodgerblue') +
    theme_classic() + 
    geom_line(data = data_frame(size = 1:1000, prop = act_prop(size)) %>% 
                  mutate(prop = max(frag_sizes$prop) * prop / max(prop)),
              aes(size), size = 0.75) +
    ylab('Proportion of total fragments') +
    xlab('Fragment size (bp)')


frag_sizes %>% 
    filter(type == 'predicted') %>% 
    ggplot(aes(y = prop)) + 
    geom_bar(aes(size - 25), stat = 'identity', position = 'dodge', 
             fill = 'red') +
    theme_classic() + 
    geom_line(data = data_frame(size = 1:1000, prop = pred_prop(size)) %>% 
                  mutate(prop = max(frag_sizes$prop[frag_sizes$type == 'predicted']) *
                             prop / max(prop)),
              aes(size), size = 0.75) +
    ylab('Proportion of total fragments') +
    xlab('Fragment size (bp)')



# Comparing distributions:
data_frame(x = 1:1000) %>% 
    ggplot(aes(x)) +
    theme_classic() +
    geom_line(aes(y = act_prop(x)), color = 'dodgerblue', size = 0.75) +
    geom_line(aes(y = pred_prop(x)), color = 'red', size = 0.75) +
    xlab('Fragment length (bp)') +
    ylab('Density')



sim_frags <- data_frame(frag_len = rnbinom(3e4L, size = coef(predicted_fit)[['size']], 
                                           mu = coef(predicted_fit)[['mu']])) %>% 
    mutate(
        frag_len = ifelse(frag_len == 0, 1, frag_len),
        prop = bin_prop(frag_len),
        id = paste0('f_', 1:length(frag_len))) %>% 
    group_by(prop) %>% 
    mutate(keep = rbinom(n(), 1, p = prop)) %>% 
    ungroup


sim_frags %>% 
    filter(keep == 1) %T>% 
    {N <<- nrow(.); .} %>% 
    ggplot(aes(frag_len)) +
    theme_classic() +
    geom_histogram(aes(y = ..density..), binwidth = 50, fill = 'dodgerblue') +
    annotate('text', x = 400, y = 0.002, 
             label = sprintf('paste(italic(p) == %.4f)', N / nrow(sim_frags)), 
             parse = TRUE) +
    coord_cartesian(xlim = c(0, 1000))



# sim_frags %>% 
#     ggplot(aes(factor(keep, levels = 0:1, labels = c('No', 'Yes')), frag_len)) +
#     geom_point(position = position_jitter(0.2, 0), alpha = 0.2) +
#     theme_classic() +
#     scale_y_continuous('Fragment length (bp)') +
#     scale_x_discrete('Kept in simulation?')




# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================

#           Reading in fragments

# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================


chosen_res <- c('ApeKI', 'BstBI', 'NruI-HF')

dig_frags <- lapply(setNames(chosen_res, chosen_res), function(enz) {
    fasta <- readFasta(sprintf('frags_%s.fa.gz', enz))
    sread(fasta)
})

dig_frag_sizes <- lapply(dig_frags, width)

dig_frag_df <- lapply(chosen_res, 
                      function(re){
                          data_frame(frag_len = dig_frag_sizes[[re]]) %>% 
                              mutate(prop = act_prop(frag_len), enzyme = re)
                      }) %>% 
    bind_rows

do.call(dnbinom, as.list(c(coef(actual_fit), x = 10)))
dnbinom(10, size = coef(actual_fit)[['size']], 
        prob = (coef(actual_fit)[['size']] / {coef(actual_fit)[['mu']] + coef(actual_fit)[['size']]}))
dnbinom(10, size = coef(actual_fit)[['size']], 
        prob = (coef(actual_fit)[['size']] / {coef(actual_fit)[['mu']] + coef(actual_fit)[['size']]}))
# dlnorm(1e4, meanlog = 4.7235121, sdlog = 0.6949672)



differ_curves <- function(par_name, par_vec, p_cols = c('#e41a1c','#377eb8','#4daf4a',
                                                        '#984ea3','#ff7f00')) {
    if(length(par_vec) != length(p_cols)) stop('par_vec must be same length as p_cols')
    curve(dnbinom(x, size = coef(actual_fit)[['size']], mu = coef(actual_fit)[['mu']]), 
          1, 1000, n = 1000, ylab = 'density', main = par_name)
    if (par_name == 'mu') {
        for (z in 1:length(par_vec)){
            curve(dnbinom(x, size = coef(actual_fit)[['size']], 
                          mu = coef(actual_fit)[['mu']] + par_vec[z]),
                  1, 1000, n = 1000, add = TRUE, col = p_cols[z])
        }
    } else if (par_name == 'size') {
        for (z in 1:length(par_vec)){
            curve(dnbinom(x, size = coef(actual_fit)[['size']] + par_vec[z], 
                          mu = coef(actual_fit)[['mu']]),
                  1, 1000, n = 1000, add = TRUE, col = p_cols[z])
        }
    } else {
        for (z in 1:length(par_vec)){
            curve(dnbinom(x, size = coef(actual_fit)[['size']], 
                          prob = (coef(actual_fit)[['size']] / {coef(actual_fit)[['mu']] +
                                  coef(actual_fit)[['size']]}) + par_vec[z]),
                  1, 1000, n = 1000, add = TRUE, col = p_cols[z])
        }
    }
}

par(mfrow = c(3, 1))
differ_curves('size', 1:5*1)
differ_curves('mu', 1:5*20)
differ_curves('prob', -1:-5*0.001)
par(mfrow = c(1,1))


enz <- chosen_res[2]
dig_frag_df %>% 
    filter(enzyme == enz) %>% 
    sample_frac(seq_p, weight = prop) %T>%
    {N <<- nrow(.); .} %>%
    ggplot(aes(frag_len)) +
    geom_histogram(aes(y = ..density..), binwidth = 50, fill = 'dodgerblue') +
    # annotate('text', x = 400, y = 0.002,
    #          label = sprintf('paste(italic(p) == %.4f)', 
    #                          N / nrow(filter(dig_frag_df, enzyme == enz))),
    #          parse = TRUE) +
    # coord_cartesian(xlim = c(0, 1000)) +
    theme_classic()




.test <- filter(dig_frag_df, enzyme == 'ApeKI')

filt_frags <- function(props, multiplier = 6e2) {
    .keep <- sapply(props * multiplier, function(p) rbinom(1, 1, prob = min(c(1, p))))
    return(.keep == 1)
}

# Testing multipliers:
mean(filt_frags(.test$prop, 6e2))




# Visualizing output:
data_frame(size = .test$frag_len[filt_frags(.test$prop)]) %>% 
    ggplot() +
    theme_classic() +
    geom_histogram(aes(size, ..density..), binwidth = 50, fill = 'dodgerblue', 
                   closed = 'left', center = 25) +
    geom_line(data = data_frame(size = 1:750, prop = act_prop(size)) %>% 
                  mutate(prop = 0.0055 * prop / max(prop)),
              aes(size, prop), size = 0.75)





rm(.test)

