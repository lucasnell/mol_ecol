Digest genome
================
Lucas Nell
19 March 2017

-   [Read genome](#read-genome)
-   [Make restriction enzyme data frame](#make-restriction-enzyme-data-frame)
-   [Digest genome](#digest-genome)
    -   [Accounting for missing data](#accounting-for-missing-data)
    -   [Digestion summary](#digestion-summary)
-   [Choosing enzymes and visualizing fragment sizes](#choosing-enzymes-and-visualizing-fragment-sizes)
-   [Session info and package versions](#session-info-and-package-versions)

*Updated 06 April 2017*

In this script I perform in silico digestions of the aphid genome using multiple restriction enzymes. The goal here is to figure out which enzymes to use for simulations. See [`../wr_files/digest.R`](../wr_files/digest.R) for working R objects created from this script that are used for downstream processes.

**Loading packages:**

``` r
suppressPackageStartupMessages({
    library(SimRAD)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(ggplot2)
})
```

*Note*: Installing `SimRAD` requires the following code:

``` r
source('https://bioconductor.org/biocLite.R')
biocLite('Biostrings')
biocLite('ShortRead')
biocLite('zlibbioc')
install.packages('SimRAD')
```

Read genome
===========

This converts the compressed fasta file of the aphid genome to a single string containing a randomly chosen 10% of the sequences in the file. I'm using only 10% for testing because using all sequences takes a long time and uses a lot of memory.

(See the [`README.md`](../README.md) file for why I'm including `../genome_data/` in file paths.)

``` r
set.seed(63)
genome_seq <- ref.DNAseq('../genome_data/aphid_genome.fa.gz', prop.contigs = 0.1)
```

If you're more patient than me and want to test this script on the entire genome, you can run the following code instead:

``` r
genome_seq <- ref.DNAseq('../genome_data/aphid_genome.fa.gz', subselect.contigs = FALSE)
```

Make restriction enzyme data frame
==================================

Here is a table of the restriction enzymes considered and their restriction-site sequences:

| enzyme    |   site1   |  site2 |
|:----------|:---------:|:------:|
| *ApeKI*   |   G/CAGC  | G/CTGC |
| *SbfI*    | CCTGCA/GG |        |
| *PstI*    |  CTGCA/G  |        |
| *EcoT22I* |  ATGCA/T  |        |
| *BstBI*   |  TT/CGAA  |        |
| *AscI*    | GG/CGCGCC |        |
| *BspEI*   |  T/CCGGA  |        |
| *AclI*    |  AA/CGTT  |        |
| *FspI*    |  TGC/GCA  |        |
| *MluI-HF* |  A/CGCGT  |        |
| *NruI-HF* |  TCG/CGA  |        |

Below is a data frame of the above table:

``` r
enz_df <- data_frame(enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
                               'BspEI', 'AclI', 'FspI', 'MluI-HF', 'NruI-HF'),
                    sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                                 c('CTGCA', 'G'), c('ATGCA', 'T'), c('TT', 'CGAA'), 
                                 c('GG', 'CGCGCC'), c('T', 'CCGGA'), c('AA', 'CGTT'), 
                                 c('TGC', 'GCA'), c('A', 'CGCGT'), c('TCG','CGA'))) %>% 
    arrange(enzyme)
```

Digest genome
=============

This function runs an in silico digestion of the genome, given an enzyme's site sequences as a character vector:

``` r
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
```

Running that on all digestion enzymes in `enz_df`:

``` r
enz_df <- enz_df %>% 
    mutate(digest = lapply(sites, digest_genome))
```

### Accounting for missing data

> From six sequencing lanes, we identified 809,651 sequence tags (at least five times) from one or both flanks of 654,998 of the 2.1 million *ApeKI* cut sites lying within the single copy genomic fraction.

([Elshire et al. 2011](http://dx.plos.org/10.1371/journal.pone.0019379), p 5)

From above, I'm creating an object storing the proportion of cut sites that I'll assume get sequenced:

``` r
seq_p <- 654998 / 2.1e6
```

### Digestion summary

Printing summary of each digestion, where all numbers assume `seq_p` proportion of sites get sequenced:

    ## ---    AclI    ----
    ## Loci per Mbp = 93.15 
    ## Total loci = 57,892 
    ## 
    ## ---    ApeKI    ----
    ## Loci per Mbp = 180.38 
    ## Total loci = 112,104 
    ## 
    ## ---    AscI    ----
    ## Loci per Mbp = 0.62 
    ## Total loci = 384 
    ## 
    ## ---    BspEI    ----
    ## Loci per Mbp = 19.57 
    ## Total loci = 12,161 
    ## 
    ## ---    BstBI    ----
    ## Loci per Mbp = 77.32 
    ## Total loci = 48,052 
    ## 
    ## ---    EcoT22I    ----
    ## Loci per Mbp = 105.00 
    ## Total loci = 65,257 
    ## 
    ## ---    FspI    ----
    ## Loci per Mbp = 27.01 
    ## Total loci = 16,787 
    ## 
    ## ---    MluI-HF    ----
    ## Loci per Mbp = 47.31 
    ## Total loci = 29,403 
    ## 
    ## ---    NruI-HF    ----
    ## Loci per Mbp = 23.35 
    ## Total loci = 14,510 
    ## 
    ## ---    PstI    ----
    ## Loci per Mbp = 59.82 
    ## Total loci = 37,176 
    ## 
    ## ---    SbfI    ----
    ## Loci per Mbp = 1.84 
    ## Total loci = 1,142

Choosing enzymes and visualizing fragment sizes
===============================================

From the summary above (and because all the following enzymes have overhangs), I'll use *ApeKI* as a common restriction enzyme, *BstBI* as intermediate, and *BspEI* as rare.

``` r
chosen_enz <- c('ApeKI', 'BstBI', 'BspEI')
```

Below are histograms of the fragment sizes for the genome digested with each enzyme.

![](digest_files/figure-markdown_github/plot_frag_sizes-1.png)

Session info and package versions
=================================

    ## Session info --------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.3.3 (2017-03-06)
    ##  system   x86_64, darwin13.4.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-04-06

    ## Packages ------------------------------------------------------------------

    ##  package              * version  date       source        
    ##  assertthat             0.1      2013-12-06 CRAN (R 3.3.0)
    ##  backports              1.0.5    2017-01-18 CRAN (R 3.3.2)
    ##  Biobase              * 2.34.0   2016-10-18 Bioconductor  
    ##  BiocGenerics         * 0.20.0   2016-10-18 Bioconductor  
    ##  BiocParallel         * 1.8.1    2016-10-30 Bioconductor  
    ##  Biostrings           * 2.42.1   2016-12-01 Bioconductor  
    ##  bitops                 1.0-6    2013-08-17 CRAN (R 3.3.0)
    ##  colorspace             1.3-2    2016-12-14 CRAN (R 3.3.2)
    ##  DBI                    0.6      2017-03-09 CRAN (R 3.3.2)
    ##  devtools               1.12.0   2016-06-24 CRAN (R 3.3.0)
    ##  digest                 0.6.12   2017-01-27 CRAN (R 3.3.2)
    ##  dplyr                * 0.5.0    2016-06-24 CRAN (R 3.3.0)
    ##  evaluate               0.10     2016-10-11 CRAN (R 3.3.0)
    ##  GenomeInfoDb         * 1.10.3   2017-02-07 Bioconductor  
    ##  GenomicAlignments    * 1.10.1   2017-03-18 Bioconductor  
    ##  GenomicRanges        * 1.26.4   2017-03-18 Bioconductor  
    ##  ggplot2              * 2.2.1    2016-12-30 CRAN (R 3.3.2)
    ##  gtable                 0.2.0    2016-02-26 CRAN (R 3.3.0)
    ##  highr                  0.6      2016-05-09 CRAN (R 3.3.0)
    ##  htmltools              0.3.5    2016-03-21 CRAN (R 3.3.0)
    ##  hwriter                1.3.2    2014-09-10 CRAN (R 3.3.0)
    ##  IRanges              * 2.8.2    2017-03-18 Bioconductor  
    ##  knitr                  1.15.1   2016-11-22 CRAN (R 3.3.2)
    ##  labeling               0.3      2014-08-23 CRAN (R 3.3.0)
    ##  lattice                0.20-35  2017-03-25 CRAN (R 3.3.2)
    ##  latticeExtra           0.6-28   2016-02-09 CRAN (R 3.3.0)
    ##  lazyeval               0.2.0    2016-06-12 CRAN (R 3.3.0)
    ##  magrittr               1.5      2014-11-22 CRAN (R 3.3.0)
    ##  Matrix                 1.2-8    2017-01-20 CRAN (R 3.3.3)
    ##  memoise                1.0.0    2016-01-29 CRAN (R 3.3.0)
    ##  munsell                0.4.3    2016-02-13 CRAN (R 3.3.0)
    ##  plyr                   1.8.4    2016-06-08 CRAN (R 3.3.0)
    ##  purrr                * 0.2.2    2016-06-18 CRAN (R 3.3.0)
    ##  R6                     2.2.0    2016-10-05 CRAN (R 3.3.0)
    ##  RColorBrewer           1.1-2    2014-12-07 CRAN (R 3.3.0)
    ##  Rcpp                   0.12.10  2017-03-19 CRAN (R 3.3.2)
    ##  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.3.0)
    ##  reshape2               1.4.2    2016-10-22 CRAN (R 3.3.0)
    ##  rmarkdown              1.4      2017-03-24 CRAN (R 3.3.2)
    ##  rprojroot              1.2      2017-01-16 CRAN (R 3.3.2)
    ##  Rsamtools            * 1.26.1   2016-10-22 Bioconductor  
    ##  S4Vectors            * 0.12.2   2017-03-18 Bioconductor  
    ##  scales                 0.4.1    2016-11-09 CRAN (R 3.3.2)
    ##  ShortRead            * 1.32.1   2017-03-18 Bioconductor  
    ##  SimRAD               * 0.96     2016-01-06 CRAN (R 3.3.0)
    ##  stringi                1.1.3    2017-03-21 CRAN (R 3.3.2)
    ##  stringr                1.2.0    2017-02-18 CRAN (R 3.3.2)
    ##  SummarizedExperiment * 1.4.0    2016-10-18 Bioconductor  
    ##  tibble                 1.2      2016-08-26 CRAN (R 3.3.0)
    ##  tidyr                * 0.6.1    2017-01-10 CRAN (R 3.3.2)
    ##  withr                  1.0.2    2016-06-20 CRAN (R 3.3.0)
    ##  XVector              * 0.14.1   2017-03-18 Bioconductor  
    ##  yaml                   2.1.14   2016-11-12 CRAN (R 3.3.2)
    ##  zlibbioc             * 1.20.0   2016-10-18 Bioconductor
