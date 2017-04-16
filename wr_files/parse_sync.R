#!/usr/bin/env Rscript

library(Rcpp)
library(stringr)
library(dplyr)
library(readr)
library(purrr)

sourceCpp('./wr_files/parse_sync.cpp')


cpp_split_int(line, 2)
(cpp_split_delim(line) %>% length) - 4
str_split(line, '\t')[[1]]

sf <- '/Volumes/750gb/allele_freq/apeki_small.sync'


r_parse <- function(fn) {
    file_list <- read_lines(fn) %>% 
        str_split('\t')
    
    mat12 <- file_list %>% 
        map(~ matrix(.x[1:2], nrow = 1)) %>% 
        do.call(what = rbind, args = .)
    
    dep_mat <- file_list %>% 
        map(~ (.x %>% tail(-3) %>% str_split(':') %>% 
                   do.call(what = cbind, args = .))[1:4,] %>% 
                as.integer) %>%
                # as.character) %>%
        do.call(what = rbind, args = .)
    
    out_list <- list(scaffolds = mat12[,1], 
                     # positions = mat12[,2], 
                     positions = as.integer(mat12[,2]), 
                     depths = dep_mat)
}
# file_list <- read_lines(sf) %>% 
#     str_split('\t')
# 
# out_df <- file_list %>% 
#     map_df(~ as_data_frame(matrix(.x[1:2], nrow = 1))) %>% 
#     rename_(.dots = setNames(c('V1', 'V2'), c('scaff', 'pos'))) %>% 
#     mutate(pos = as.integer(pos))
# 
# dep_df <- file_list %>% 
#     map(~ (.x %>% tail(-3) %>% str_split(':') %>% 
#                do.call(what = cbind, args = .))[1:4,] %>% 
#             as.integer) %>% 
#     do.call(what = rbind, args = .) %>% 
#     as_data_frame
# 
# 
# colnames(dep_df) <- c(
#     paste0(c('A', 'T', 'C', 'G'), 'p'),
#     lapply(1:((ncol(dep_df) / 4) - 1), 
#            function(i) paste0(c('A', 'T', 'C', 'G'), i)) %>% 
#         c(recursive = TRUE))
# 
# out_df <- bind_cols(out_df, dep_df)


# Total # lines in final file: 9,795,971
# In this smaller file: 90,000
system.time(r <- r_parse(sf))
system.time(c <- cpp_parse_sync(sf, 10))

identical(r, c)
# str(r); str(c)



