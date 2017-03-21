#' ---
#' title: 'Digest genome'
#' author: 'Lucas Nell'
#' date: '`r format(Sys.Date())`'
#' output: github_document
#' ---
#' 
#' 
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


#+ packages, echo = FALSE
suppressPackageStartupMessages(
    suppressWarnings({
        library(SimRAD)
        library(tidyverse)
        library(fitdistrplus)
    })
)


#+ 
set.seed(699)
genome_seq <- ref.DNAseq('aphid_genome.fa.gz', prop.contigs = 0.1)
paste(substr(genome_seq, 1, 10), '...', 
      substr(genome_seq, nchar(genome_seq) - 9, nchar(genome_seq)))




re_df <- data_frame(enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI'),
                    sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                                 c('CTGCA', 'G'), c('ATGCA', 'T'), 
                                 c('TT', 'CGAA')))


digest_enzyme <- function(enzyme_sites, dna_seq = genome_seq) {
    if (is.list(enzyme_sites)) {
        enzyme_sites <- enzyme_sites[[1]]
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


re_df <- re_df %>% 
    mutate(digest = lapply(sites, digest_enzyme))



# # Restriction enzyme: ApeKI
# apeki_dig <- digest_enzyme(c('G', 'CAGC', 'G', 'CTGC'))
# 
# # SbfI
# sbfi_dig <- digest_enzyme(c('CCTGCA', 'GG'))
# 
# # PstI 
# psti_dig <- digest_enzyme(c('CTGCA', 'G'))
# 
# # EcoT22I 
# ecot22i_dig <- digest_enzyme(c('ATGCA', 'T'))

enz <- 'BstBI'
dig <- digest_enzyme(c('TT', 'CGAA'))
cat('Total loci for', enz, '=', 
    format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))

enz <- 'AscI'
dig <- digest_enzyme(c('GG', 'CGCGCC'))
cat('Total loci for', enz, '=', 
    format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))

enz <- 'BspEI'
dig <- digest_enzyme(c('T', 'CCGGA'))
cat('Total loci for', enz, '=', 
    format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))

enz <- 'AclI'
dig <- digest_enzyme(c('AA', 'CGTT'))
cat('Total loci for', enz, '=', 
    format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))

# enz <- 'FspI'
# dig <- digest_enzyme(c('TGC', 'GCA'))
# cat('Total loci for', enz, '=', 
#     format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))
# # Not considered because "Ligation is 25% -75%."

enz <- 'MluI-HF'
dig <- digest_enzyme(c('A', 'CGCGT'))
cat('Total loci for', enz, '=', 
    format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))

enz <- 'NruI-HF'
dig <- digest_enzyme(c('TCG','CGA'))
cat('Total loci for', enz, '=', 
    format(seq_p * length(dig) * 10, digits = 0, big.mark = ',', scientific = FALSE))


#' 
#' > From six sequencing lanes, we identified 809,651 sequence tags (at least five times)
#' > from one or both flanks of 654,998 of the 2.1 million ApeKI cut sites lying within
#' > the single copy genomic fraction.
#' 
#' (Elshire et al. 2011, p 5)
#' 
#' Below is the proportion of cut sites that get sequenced:
#' 

seq_p <- 654998 / 2.1e6


#'
#'  Printing summary of each digestion:
#'  
#+ dig_summary, echo = FALSE
z <- apply(re_df, 1, 
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
           }); rm(z)





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




