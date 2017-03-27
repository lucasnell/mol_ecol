#' ---
#' title: "Choosing the sequencing simulator"
#' author: "Lucas Nell"
#' date: "19 March 2017"
#' output: 
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#'     fig_width: 8
#'     fig_height: 6
#' ---
#' 
#' *Updated `r format(Sys.Date(), '%d %B %Y')`*
#' 
#' 
#' Tables are from the following paper (link 
#' [here](http://www.nature.com/nrg/journal/v17/n8/full/nrg.2016.57.html)):
#' 
#' Escalona, M., S. Rocha, and D. Posada. 2016. A comparison of tools for the simulation 
#' of genomic next-generation sequencing data. *Nature Reviews Genetics* __17__:459–469.
#' 
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages({
    library(dplyr)
    library(readxl)
})

#' 
#' ## Reading Excel sheets
#' 

sheet_names <- c(paste0('table', 2:4), paste0('tab', 2:4, '_meta'))

sheet_list <- lapply(sheet_names, function(sh) read_excel('ngs_simulators.xlsx', sh))

meta_data <- bind_rows(sheet_list[4:6]) %>% 
    arrange(abbrev) %>% 
    as.data.frame

# Joining all tables together (not including metadata):
all_tables <- left_join(sheet_list[[1]], sheet_list[[2]], by = 'Simulators') %>% 
    left_join(., sheet_list[[3]], by = 'Simulators') %>% 
    arrange(Simulators) %>% 
    # Removing spaces from column names:
    rename_(.dots = setNames(sprintf('`%s`', colnames(.)), 
                             gsub(' ', '_', colnames(.)))) %>%
    # Keeping only columns I'll find useful:
    select_(
        .dots = c(
            'Simulators',
            # From Table 2:
            'Technology', 'G_vs_M', 'Run_types', 'PCR', 'QS', 'out_RE', 'AL', 'FO', 
            # From Table 3:
            'Programming_language', 'Operating_system', 'Processing', 'Open_source', 
            # From Table 4:
            'SNPs'
        )
    )




meta_data

all_tables %>% 
    filter(
        # grepl('Illumina', Technology),
        # grepl('SE', Run_types),
        # PCR == 'Yes',
        grepl('Mac', Operating_system),
        # Open_source == 'Yes',
        SNPs == 'Yes'
        ) %>% 
    # select(-Technology, -Run_types, -PCR, -Operating_system, -Open_source, -SNPs) %>%
    # select(Simulators) %>%
    as.data.frame

# For filtering by just one at a time:
eval_list <- list(
    illumina = expression(grepl('Illumina', Technology)),
    single_end = expression(grepl('SE', Run_types)),
    pcr = expression(PCR == 'Yes'),
    mac = expression(grepl('Mac', Operating_system)),
    open = expression(Open_source == 'Yes'),
    snps = expression(SNPs == 'Yes')
)

with(all_tables, filter(all_tables, eval(eval_list[['pcr']]))) %>% as.data.frame

with(all_tables, filter(all_tables, eval(eval_list[['snps']]))) %>% as.data.frame








