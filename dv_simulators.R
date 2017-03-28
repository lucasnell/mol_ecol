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
#' 
#' ## Reading Excel sheets
#' 
sheet_names <- c(paste0('table', 2:4), paste0('tab', 2:4, '_meta'))
sheet_list <- lapply(sheet_names, 
                     function(sh) read_excel('./bg_data/ngs_simulators.xlsx', sh))
#' 
#' Joining together metadata tables:
#' 
meta_data <- bind_rows(sheet_list[4:6]) %>% 
    arrange(abbrev) %>% 
    as.data.frame

#' 
#' Joining all other tables together:
#' 
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
#' 
#' 
#' Abbreviations for some of the `all_tables` column names:
#' 
meta_data %>% filter(abbrev %in% c(colnames(all_tables), 'RE'))
#' 
#' 
#' 
#' # Choosing a simulator
#' 
#' 
#' I need it to do the following:
#' 
#' * Work on a Mac
#' * Simulate single-end, Illumina reads
#' * Simulate the PCR step
#' * Be open source
#' 
#' I would also like it to simulate SNPs from a reference genome, so I'll do this 
#' filtering here:
#' 
all_tables %>% 
    filter(
        grepl('Illumina', Technology),
        grepl('SE', Run_types),
        PCR == 'Yes',
        grepl('Mac', Operating_system),
        Open_source == 'Yes',
        SNPs == 'Yes'
        ) %>% 
    select(Simulators)
#' 
#' 
#' When I tried out Grinder, it was not very intuitive, especially the SNPs part, so I
#' decided to remove that filter:
#' 
all_tables %>% 
    filter(
        grepl('Illumina', Technology),
        grepl('SE', Run_types),
        PCR == 'Yes',
        grepl('Mac', Operating_system),
        Open_source == 'Yes'
    ) %>% 
    select(Simulators)
#' 
#' 
#' This includes the ART simulator. I will use this one for read simulation, and 
#' generate SNPs myself.
#' 
#' 