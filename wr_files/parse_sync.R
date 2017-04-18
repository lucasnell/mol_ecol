

# Parsing sync files output from PoPoolation2


library(Rcpp)
library(stringr)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))


sourceCpp('./wr_files/parse_sync.cpp')



genome <- cpp_get_genome('./genome_data/aphid_genome.fa') %>% 
    do.call(what = cbind, args = .) %>% 
    as_data_frame %>% 
    rename(name = names, size = sizes) %>% 
    mutate(size = as.integer(size)) %>% 
    arrange(desc(size))


# # Took ~300 seconds
# system.time(sync <- cpp_parse_sync_all_info('/Volumes/750gb/allele_freq/apeki.sync',
#                                             n_samps = 10))
# rm(sync); invisible(gc())






parse_sync <- function(fn, n_samps = 10, filter_same = TRUE) {
    
    sync_list <- cpp_parse_sync_concise(fn, n_samps)

    pos_df <- data_frame(scaf = sync_list$scaffolds) %>% 
        mutate(
            pos = sync_list$positions,
            dep = rowSums(sync_list$pool_depths)
        )
    
    allele_df <- as_data_frame(sync_list$samp_alleles)
    colnames(allele_df) <- paste0('s_', 
                                  sprintf(1:n_samps, 
                                          fmt = paste0('%0', nchar(n_samps), 'i')))
    out_df <- bind_cols(pos_df, allele_df)
    
    if (filter_same) {
        out_df <- filter(out_df, !sync_list$all_equal)
    }
    
    out_df <- gather(out_df, 'samp', 'freq', starts_with("s_"))
    
    return(out_df)
}



# Figure of read depths (commented bc it takes ~5-6 minutes)
system.time(sync <- parse_sync('/Volumes/750gb/allele_freq/apeki.sync', 10, FALSE))

ss <- sync %>%
    filter(samp == 's_01', dep > 0)
dep_plot <- ss %>%
    ggplot(aes(dep, ..density..)) +
    geom_histogram(bins = 100, fill = 'dodgerblue') +
    xlab('Read depth') +
    ylab('Density')

ggsave('./_paper/figs/dep_p.pdf', dep_plot, width = 6, height = 2.5)



# # Takes ~5-6 minutes
# system.time(sync <- parse_sync('/Volumes/750gb/allele_freq/apeki.sync'))
# save(sync, file = 'sync.RData')

load('sync.RData')

sync <- sync %>%
    mutate_if(is.character, factor)

sync %>% filter(samp == 's_01')

sync %>% 
    group_by(samp) %>% 
    summarize(p = mean(freq)) %>% 
    mutate(p = round(p / sum(p), 4)) %>% 
    knitr::kable(format = 'latex')



prop_p <- sync %>% 
    group_by(samp) %>% 
    summarize(p = mean(freq)) %>% 
    mutate(p = p / sum(p), samp = factor(as.integer(samp))) %>% 
    ggplot(aes(samp, p)) +
    geom_bar(stat = 'identity', fill = 'dodgerblue') +
    geom_hline(yintercept = 0.1, linetype = 3) +
    ylab('Estimated proportion') +
    xlab('Sample')
prop_p
ggsave('./_paper/figs/prop_p.pdf', prop_p, width = 6, height = 3)




data_frame(X = seq(0, 75, 5)) %>% 
    mutate(sd_X = sapply(X, function(x) {
        m <- (sync %>% 
                  filter(dep > x) %>%
                  group_by(samp) %>% 
                  summarize(m = mean(freq)) %>% 
                  mutate(m = m / sum(m)))$m
        return(sd(m))
    })) %>% 
    ggplot(aes(X, sd_X)) +
    geom_line() +
    ylab('SD of proportion estimates') +
    xlab('Depth cutoff (bp)')



sync %>% 
    # filter(scaf == 'GL349630.1') %>% 
    group_by(pos) %>% 
    mutate(freq = freq / sum(freq)) %>% 
    ungroup %>% 
    ggplot(aes(samp, freq, fill = samp)) + 
    geom_violin(scale = 'count') +
    # geom_point(alpha = 0.1, position = position_jitter(0.25, 0), size = 0.5) +
    # geom_line() +
    # facet_wrap(~ samp, nrow = 2)
    theme(legend.position = 'none') +
    ylab('Frequency') +
    xlab('Sample')














# ---------------------------------------
# ---------------------------------------

# Plotting depths

# ---------------------------------------
# ---------------------------------------
# 
# expand_scaf <- function(scaf_df, genome_df = genome) {
#     options(stringsAsFactors = FALSE)
#     this_scaf <- scaf_df$scaf[1]
#     max_pos <- genome_df$size[genome_df$name == this_scaf]
#     out_df <- scaf_df
#     if (min(scaf_df$pos) > 1) {
#         min_row <- expand.grid(scaf = this_scaf, pos = 1:(min(scaf_df$pos) - 1),
#                                dep = 0, samp = unique(scaf_df$samp), freq = 0) %>%
#             mutate_at(vars(scaf, samp), funs(as.character))
#         out_df <- bind_rows(min_row, out_df)
#     }
#     if (max(scaf_df$pos) < max_pos) {
#         max_row <- expand.grid(scaf = this_scaf, pos = (max(scaf_df$pos) + 1):max_pos,
#                                dep = 0, samp = unique(scaf_df$samp), freq = 0) %>%
#             mutate_at(vars(scaf, samp), funs(as.character))
#         out_df <- bind_rows(out_df, max_row)
#     }
#     out_df <- as_data_frame(out_df)
#     return(out_df)
# }
# 
# exp_sync <- sync %>% 
#     filter(scaf %in% genome$names[1:10]) %>%
#     split(.$scaf) %>% 
#     map(expand_scaf) %>% 
#     bind_rows
# 
# plot(1:nrow(genome), genome$sizes)
# 
# library(RcppRoll)
# # Read depth
# exp_sync %>% 
#     filter(samp == 's_01') %>% 
#     group_by(scaf) %>% 
#     summarize(dep = list(roll_mean(dep, n = 1e3L, by = 500L, align = 'left')),
#               g_pos = list(1:length(dep[[1]]))) %>% 
#     unnest %>% 
#     mutate(scaf = factor(scaf)) %>% 
#     ggplot(aes(g_pos / 1e3, dep, color = scaf)) +
#     # geom_point(alpha = 0.1)
#     geom_line() +
#     facet_wrap(~ scaf, nrow = 1, scales = 'free_x')












