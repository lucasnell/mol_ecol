
# Plot of sampling error for samples taken from population of binary genotypes

library(dplyr)
library(ggplot2)

pop <- c(rep(0, 500), rep(1, 500))

samp_ns <- c(25, 50, 100)
set.seed(9)
samp_df <- lapply(samp_ns, 
                  function(x) sapply(1:1000, function(i) mean(sample(pop, x)))) %>% 
    do.call(what = cbind, args = .) %>% 
    as_data_frame %>% 
    rename_(.dots = setNames(paste0('V', 1:3), samp_ns)) %>% 
    tidyr::gather('N', mu)

sampling_p <- samp_df %>% 
    filter(N != '5') %>% 
    mutate(N = factor(N, levels = samp_ns)) %>% 
    ggplot(aes(mu, fill = N)) +
    geom_density(color = NA, alpha = 0.7) +
    scale_fill_manual(values = c('#d7191c','#fdae61','#2c7bb6')) +
    xlab(expression(bar(x) ~ estimates)) +
    ylab('Density') +
    theme(legend.position = c(0.1, 0.8))
sampling_p
ggsave('./_paper/figs/sampling_p.pdf', sampling_p, width = 6, height = 4)

