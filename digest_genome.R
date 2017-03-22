#' ---
#' title: 'Digest genome'
#' author: 'Lucas Nell'
#' date: '19 March 2017'
#' output: github_document
#' ---
#' 
#' *Updated `r format(Sys.Date(), '%d %B %Y')`*
#' 
#' 
#' In this script I perform in silico digestions of the aphid genome using multiple 
#' restriction enzymes.
#' Once enzymes are chosen, I write the resulting fragments to fasta files.
#' 
#' 
#' __Loading packages:__
#' 
#+ packages
suppressPackageStartupMessages(
    suppressWarnings({
        library(SimRAD)
        library(dplyr)
        library(purrr)
        library(tidyr)
        library(ggplot2)
    })
)
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
#' containing the sequences in the file. (It should take ~20 seconds.)
#' 
#+ make_genome
genome_seq <- ref.DNAseq('aphid_genome.fa.gz', subselect.contigs = FALSE)
#' 
#' 
#' If you want to test this script without using as much RAM or time, you can run the
#' following code instead:
#' 
#' ```{r, eval = FALSE}
#' set.seed(1)
#' genome_seq <- ref.DNAseq('aphid_genome.fa.gz', prop.contigs = 0.1)
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
#+ make_re_df
re_df <- data_frame(enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
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
digest_enzyme <- function(enzyme_sites, dna_seq = genome_seq) {
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
#' Running that on all digestion enzymes in `re_df` (__Warning:__ this takes ~ 4.5 
#' minutes and can use > 4GB RAM):
#' 
#+ run_digest
re_df <- re_df %>% 
    mutate(digest = lapply(sites, digest_enzyme))

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
#' ### Digestion summary
#' 
#' Printing summary of each digestion, where all numbers assume `seq_p` proportion of
#' sites get sequenced:
#'  
#+ dig_summary, echo = FALSE
z <- apply(re_df, 1, 
           function(i){
               dig_i <- i$digest
               name_i <- i$enzyme
               loci_density_i <- seq_p * length(dig_i) / (nchar(genome_seq) / 1e6)
               loci_total_i <- seq_p * length(dig_i)
               cat('---   ', name_i, '   ----\n')
               cat(sprintf('Loci per Mbp = %.2f', loci_density_i), '\n')
               cat(sprintf('Total loci = %s', 
                           format(loci_total_i, big.mark = ',', 
                                  digits = 0, scientific = FALSE)), '\n\n')
               return(NULL)
           }); rm(z)
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
plot_df <- re_df %>% 
    filter(enzyme %in% chosen_res) %>% 
    split(.$enzyme) %>% 
    map_df(~ data_frame(enzyme = .x$enzyme, frag_len = nchar(.x$digest[[1]]))) %>% 
    mutate(enzyme = factor(enzyme, levels = chosen_res))
plot_df %>% 
    ggplot(aes(frag_len / 1e3)) +
    theme_classic() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = 14, face = 'bold.italic')) +
    geom_histogram(aes(y = ..density.., fill = enzyme), bins = 100) +
    facet_grid(enzyme ~ ., scales = 'free') +
    ylab('Density') +
    scale_x_log10('Fragment length (Kbp)', breaks = c(0.1, 1, 10, 50)) +
    scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb'), guide = FALSE)
rm(plot_df); invisible(gc(verbose = FALSE))
#' 
#' 
#' 
#' 
#' 
#' # Writing to fasta files
#' 
#' The `DNAStringSet` function is from the `Biostrings` package, and `writeFasta` is
#' from the `ShortRead` package. Both of these packages should already be loaded from
#' by `SimRAD`.
#' 
#' The below code makes a list of `DNAStringSet` objects with individual sequence names 
#' set to `seq_X`, where `X` goes from 1 to the number of sequences.
#' 
#+ make_write_list
write_list <- re_df %>% 
    filter(enzyme %in% chosen_res) %>% 
    split(.$enzyme) %>% 
    map(~ DNAStringSet(.x$digest[[1]])) %>% 
    map(~ magrittr::set_names(.x, paste0('seq_', 1:length(.x))))
write_list
#' 
#' To save some RAM before writing, I'm going to dump two large objects:
#+ dump_objs
rm(re_df, genome_seq)
invisible(gc(verbose = FALSE))
#' 
#' 
#' 
#' Now I write each `DNAStringSet` object from `write_list` to a compressed fasta file 
#' (this took ~5.5 mins on my computer):
#' 
#+ write_fastas, eval = FALSE
for (enz in chosen_res) {
    writeFasta(write_list[[enz]], file = sprintf('frags_%s.fa.gz', enz), mode = 'w',
               compress = 'gzip')
    cat(sprintf('%s file finished', enz), '\n')
}
#' 
#' 
#' # Session info and package versions
#' 
#' 
#+ session_info, echo = FALSE
devtools::session_info()




