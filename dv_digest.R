#' ---
#' title: 'Digest genome'
#' author: 'Lucas Nell'
#' date: '19 March 2017'
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
#' In this script I perform in silico digestions of the aphid genome using multiple 
#' restriction enzymes.
#' The goal here is to figure out which enzymes to use for simulations.
#' See [`wr_digest.R`](./wr_digest.R) for working R objects created from this script 
#' that are used for downstream processes.
#' 
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages({
    library(SimRAD)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(ggplot2)
})
#' 
#+ set_theme, echo = FALSE
# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))
#' 
#' *Note*: Installing `SimRAD` requires the following code:
#' 
#' ```{r install_SimRAD, eval = FALSE}
#' source('https://bioconductor.org/biocLite.R')
#' biocLite('Biostrings')
#' biocLite('ShortRead')
#' biocLite('zlibbioc')
#' install.packages('SimRAD')
#' ```
#' 
#' 
#' # Read genome
#' 
#' This converts the compressed fasta file of the aphid genome to a single string
#' containing a randomly chosen 10% of the sequences in the file.
#' I'm using only 10% for testing because using all sequences takes a long time and 
#' uses a lot of memory.
#' 
#' (See the [`README.md`](./README.md) file for why I'm including `./genome_data/` in 
#' file paths.)
#' 
#+ make_genome
set.seed(63)
genome_seq <- ref.DNAseq('./genome_data/aphid_genome.fa.gz', prop.contigs = 0.1)
#' 
#' 
#' If you're more patient than me and want to test this script on the entire genome, 
#' you can run the following code instead:
#' 
#' ```{r, eval = FALSE}
#' genome_seq <- ref.DNAseq('./genome_data/aphid_genome.fa.gz', subselect.contigs = FALSE)
#' ```
#' 
#' 
#' # Make restriction enzyme data frame
#' 
#' 
#' Here is a table of the restriction enzymes considered and their restriction-site
#' sequences:
#' 
#' 
#+ display_table, echo = FALSE
disp_mat <- data_frame(enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
                                  'BspEI', 'AclI', 'FspI', 'MluI-HF', 'NruI-HF'),
                       sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                                    c('CTGCA', 'G'), c('ATGCA', 'T'), c('TT', 'CGAA'), 
                                    c('GG', 'CGCGCC'), c('T', 'CCGGA'), c('AA', 'CGTT'), 
                                    c('TGC', 'GCA'), c('A', 'CGCGT'), c('TCG','CGA'))) %>%
    mutate(sites = sapply(sites, function(s) paste0(s, collapse = ','))) %>% 
    separate(sites, paste0('s', 1:4), fill = 'right') %>% 
    unite('site1', s1:s2, sep = '/') %>% 
    unite('site2', s3:s4, sep = '/') %>% 
    mutate(site2 = ifelse(site2 == 'NA/NA', '', site2),
           enzyme = sprintf('*%s*', enzyme))
knitr::kable(disp_mat, format = 'markdown', row.names = FALSE, escape = TRUE, 
             align = c('l', 'c', 'c'))
#' 
#' 
#' 
#' Below is a data frame of the above table:
#' 
#+ make_enz_df
enz_df <- data_frame(enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
                               'BspEI', 'AclI', 'FspI', 'MluI-HF', 'NruI-HF'),
                    sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                                 c('CTGCA', 'G'), c('ATGCA', 'T'), c('TT', 'CGAA'), 
                                 c('GG', 'CGCGCC'), c('T', 'CCGGA'), c('AA', 'CGTT'), 
                                 c('TGC', 'GCA'), c('A', 'CGCGT'), c('TCG','CGA'))) %>% 
    filter(!enzyme %in% c('SbfI', 'PstI', 'EcoT22I', 'AscI', 'BspEI', 'FspI'))
#' 
#' 
#' 
#' The restriction enzymes below were filtered out (reasons and link to referencing site 
#' follow):
#' 
#' - *SbfI*: Not blocked by CpG methylase
#'   ([link](https://www.neb.com/products/r0642-sbfi))
#' - *PstI*: Not blocked by CpG methylase
#'   ([link](https://www.neb.com/products/r0140-psti))
#' - *EcoT22I*: "... not sensitive to dam, dcm, or CG methylation" 
#'   ([link](https://tools.thermofisher.com/content/sfs/manuals/15240501.pdf))
#' - *AscI*: "AscI is strongly inhibited by NaCl and ammonium acetate"
#'   ([link](https://www.neb.com/products/r0558-asci))
#' - *BspEI*: Only impaired by CpG methylase
#'   ([link](https://www.neb.com/products/r0540-bspei))
#' - *FspI*: "Ligation is 25%–75%."
#'   ([link](https://www.neb.com/products/r0135-fspi))
#' 
#' 
#' # Digest genome
#' 
#' This function runs an in silico digestion of the genome, given an enzyme's site 
#' sequences as a character vector:
#' 
#+ digest_function
digest_genome <- function(enzyme_sites, dna_seq = genome_seq) {
    if (is.list(enzyme_sites)) {
        enzyme_sites <- enzyme_sites[[1]]
    }
    if ((length(enzyme_sites) %% 2) != 0) {
        stop(paste('enzyme_sites argument must be of even length, separately providing',
                   'sequences for 5 prime and 3 prime sides of each site'))
    }
    names(enzyme_sites) <- 
        c('cut_site_5prime1', 'cut_site_3prime1', 
          'cut_site_5prime2', 'cut_site_3prime2',  
          'cut_site_5prime3', 'cut_site_3prime3', 
          'cut_site_5prime4', 'cut_site_3prime4')[1:length(enzyme_sites)]
    
    call_list <- as.list(c(DNAseq = dna_seq, verbose = FALSE, enzyme_sites))
    
    dig <- do.call(insilico.digest, call_list)
    
    return(dig)
}



#' 
#' Running that on all digestion enzymes in `enz_df`:
#' 
#+ run_digest
enz_df <- enz_df %>% 
    mutate(digest = lapply(sites, digest_genome))

#' 
#' ### Accounting for missing data
#' 
#' > From six sequencing lanes, we identified 809,651 sequence tags (at least five times)
#' > from one or both flanks of 654,998 of the 2.1 million *ApeKI* cut sites lying within
#' > the single copy genomic fraction.
#' 
#' ([Elshire et al. 2011](http://dx.plos.org/10.1371/journal.pone.0019379), p 5)
#' 
#' From above, I'm creating an object storing the proportion of cut sites that 
#' I'll assume get sequenced:
#' 
#+ make_seq_p
seq_p <- 654998 / 2.1e6
#' 
#' 
#' 
#' ### Digestion summary
#' 
#' Printing summary of each digestion, where all numbers assume `seq_p` proportion of
#' sites get sequenced:
#'  
#+ dig_summary, echo = FALSE
invisible(apply(enz_df, 1, 
                function(i){
                    dig_i <- i$digest
                    name_i <- i$enzyme
                    loci_density_i <- seq_p * length(dig_i) / (nchar(genome_seq) / 1e6)
                    loci_total_i <- seq_p * length(dig_i) * 10
                    cat('---   ', name_i, '   ----\n')
                    cat(sprintf('Loci per Mbp = %.2f', loci_density_i), '\n')
                    cat(sprintf('Total loci = %s', 
                                format(loci_total_i, big.mark = ',', 
                                       digits = 0, scientific = FALSE)), '\n\n')
                    return(NULL)
                }))
#' 
#' 
#' # Choosing enzymes and visualizing fragment sizes
#' 
#' 
#' From the summary above, I'll use *ApeKI* as a common restriction enzyme, 
#' *BstBI* as intermediate, and
#' *NruI-HF* as rare.
#' 
#+ make_chosen_res
chosen_res <- c('ApeKI', 'BstBI', 'NruI-HF')
#' 
#' Below are histograms of the fragment sizes for the genome digested with each enzyme.
#' 
#+ plot_frag_sizes, echo = FALSE
plot_df <- enz_df %>% 
    filter(enzyme %in% chosen_res) %>% 
    split(.$enzyme) %>% 
    map_df(~ data_frame(enzyme = .x$enzyme, frag_len = nchar(.x$digest[[1]]))) %>% 
    mutate(enzyme = factor(enzyme, levels = chosen_res))
plot_df %>% 
    ggplot(aes(frag_len / 1e3)) +
    theme(strip.text = element_text(size = 14, face = 'bold.italic')) +
    geom_histogram(aes(y = ..density.., fill = enzyme), bins = 100) +
    facet_grid(enzyme ~ ., scales = 'free') +
    ylab('Density') +
    scale_x_log10('Fragment length (Kbp)', breaks = c(0.1, 1, 10, 50)) +
    scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb'), guide = FALSE)
#' 
#' 
#' 
#' 
#' 
#' # Session info and package versions
#' 
#' 
#+ session_info, echo = FALSE
devtools::session_info()
