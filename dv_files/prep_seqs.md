Prepare sequences for sequencing simulation
================
Lucas Nell
28 March 2017

-   [References](#references)
-   [Function](#function)
-   [Testing function](#testing-function)
-   [Session info](#session-info)

*Updated 06 April 2017*

This script tests the function that creates strands to be used for sequencing simulation. It is to be used on digested fragments. It removes regions from fragments that will not be sequenced (those that are far from restriction enzyme cut sites). For each cut-site-adjacent area, it returns a forward and reverse strand, the latter of which is returned as its reverse complement. The reverse complement is returned because I will be simulating sequencing from fragments unidirectionally in the forward direction. I chose to do unidirectional sequencing because this more similarly simulates how both forward and reverse fragment ends can be assigned sequencing primers.

**Loading packages:**

``` r
suppressPackageStartupMessages(library(ShortRead))
```

References
==========

These are the papers I read to develop the below function:

Davey, J. W., P. A. Hohenlohe, P. D. Etter, J. Q. Boone, J. M. Catchen, and M. L. Blaxter. 2011. Genome-wide genetic marker discovery and genotyping using next-generation sequencing. *Nature Reviews Genetics* **12**:499–510. Available from [this link](http://dx.doi.org/10.1038/nrg3012).

Elshire, R. J., J. C. Glaubitz, Q. Sun, J. A. Poland, K. Kawamoto, E. S. Buckler, and S. E. Mitchell. 2011. A robust, simple genotyping-by-sequencing (GBS) approach for high diversity species. *PLOS ONE* **6**:e19379. Available from [this link](http://dx.plos.org/10.1371/journal.pone.0019379).

Function
========

The `prep_seqs` removes portions of fragments that are far from restriction enzyme cut sites. A `DNAStringSet` object of sequence fragments is the required argument. Read length (defaults to 100bp) and barcode length (defaults to 4) are optional arguments. It returns a `DNAStringSet` object with the new sequences as described above. The returned set should have twice as many sequences as the input set of sequences because `prep_seqs` outputs a sequence for each sequencing strand.

``` r
prep_seqs <- function(dna_ss, read_len = 100, bc_len = 4) {
    # Length to cut fragment to is simply the read length minus the barcode length and
    # the restriction enzyme's overhang length
    cut_len <- as.integer(read_len - bc_len)
    # Extracting character vector from input DNAStringSet
    character_seqs <- as.character(dna_ss)
    # Inner function that removes areas far from cut sites from one sequence
    .one_rm <- function(.s) {
        .seq_len <- as.integer(nchar(.s))
        if (.seq_len < cut_len) {
            f_strand <- r_strand <- .s
        } else {
            f_strand <- .Internal(substr(.s, 1L, cut_len))
            r_strand <- .Internal(substr(.s, .seq_len - cut_len + 1L, .seq_len))
        }
        return(matrix(c(f_strand, r_strand), nrow = 1))
    }
    # applying that inner function to each string in the character vector
    mat_list <- lapply(character_seqs, .one_rm)
    cut_seq_char <- do.call(rbind, mat_list)
    # Creating DNAStringSet objects for forward and reverse strands.
    # Because I am doing unidirectional (forward-only) sequencing, I need to make the 
    # reverse strands reverse complements of the forward one.
    f_strands <- DNAStringSet(cut_seq_char[,1])
    r_strands <- reverseComplement(DNAStringSet(cut_seq_char[,2]))
    # Now I combine them
    cut_dna_ss <- append(f_strands, r_strands)
    return(cut_dna_ss)
}
```

Testing function
================

I am `source`-ing `../wr_files/size_filter.R` to use those objects to first filter the fragments by size before removing faraway sequences.

``` r
source('../wr_files/size_filter.R')
```

Now I read the fasta file and do the filtering by size.

``` r
test_fasta <- sread(readFasta('../genome_data/frags_BstBI.fa.gz'))
filt_fasta <- size_filter(test_fasta)
```

Below is a run of `prep_seqs` including the run time and output.

``` r
system.time({ps_test <- prep_seqs(filt_fasta)})
```

    ##    user  system elapsed 
    ##   0.364   0.039   0.405

``` r
ps_test
```

    ##   A DNAStringSet instance of length 34906
    ##         width seq
    ##     [1]    96 CGAACTAGAGTAACGTCGTACTTGCACTT...ATCAAAAATTAACTTTTTTGGAACTCAGA
    ##     [2]    96 CGAATATTATTAAAATTTACTAACAAAAA...TTTCGTGCAATTAGATCAATAAGCGGAAA
    ##     [3]    96 CGAAACATTTTACGACCCTGAAATAGGGT...TCCTAGTTTTGGATGGATCTGTTTTAGTG
    ##     [4]     7 CGAAATT
    ##     [5]    96 CGAAATTGCCATATAATTAGAAATGAATG...GGGGTTTACGTTTAAATGTGTTATTTGTG
    ##     ...   ... ...
    ## [34902]    42 AAAATATGGTATTAGTGGTATAATTCAGCCCCCTGCCATTCG
    ## [34903]    68 AATAAATGAGAATTAAAATATTCCGTATA...ATATATTCTTAACAGTGAATATATATTCG
    ## [34904]    96 AATTAGTTACGTAATATTTGTAATGTATT...ATAGTTAAAGGAGAAACCCGTATATTTGT
    ## [34905]    12 AAAATCATTTCG
    ## [34906]    96 GGAGAGCTAGTCCTATACTATAAAATGAA...TCTCGTGGTGGAAAAAATTAAACACCGCT

Session info
============

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
    ##  backports              1.0.5    2017-01-18 CRAN (R 3.3.2)
    ##  Biobase              * 2.34.0   2016-10-18 Bioconductor  
    ##  BiocGenerics         * 0.20.0   2016-10-18 Bioconductor  
    ##  BiocParallel         * 1.8.1    2016-10-30 Bioconductor  
    ##  Biostrings           * 2.42.1   2016-12-01 Bioconductor  
    ##  bitops                 1.0-6    2013-08-17 CRAN (R 3.3.0)
    ##  devtools               1.12.0   2016-06-24 CRAN (R 3.3.0)
    ##  digest                 0.6.12   2017-01-27 CRAN (R 3.3.2)
    ##  evaluate               0.10     2016-10-11 CRAN (R 3.3.0)
    ##  GenomeInfoDb         * 1.10.3   2017-02-07 Bioconductor  
    ##  GenomicAlignments    * 1.10.1   2017-03-18 Bioconductor  
    ##  GenomicRanges        * 1.26.4   2017-03-18 Bioconductor  
    ##  htmltools              0.3.5    2016-03-21 CRAN (R 3.3.0)
    ##  hwriter                1.3.2    2014-09-10 CRAN (R 3.3.0)
    ##  IRanges              * 2.8.2    2017-03-18 Bioconductor  
    ##  knitr                  1.15.1   2016-11-22 CRAN (R 3.3.2)
    ##  lattice                0.20-35  2017-03-25 CRAN (R 3.3.2)
    ##  latticeExtra           0.6-28   2016-02-09 CRAN (R 3.3.0)
    ##  magrittr               1.5      2014-11-22 CRAN (R 3.3.0)
    ##  Matrix                 1.2-8    2017-01-20 CRAN (R 3.3.3)
    ##  memoise                1.0.0    2016-01-29 CRAN (R 3.3.0)
    ##  RColorBrewer           1.1-2    2014-12-07 CRAN (R 3.3.0)
    ##  Rcpp                   0.12.10  2017-03-19 CRAN (R 3.3.2)
    ##  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.3.0)
    ##  rmarkdown              1.4      2017-03-24 CRAN (R 3.3.2)
    ##  rprojroot              1.2      2017-01-16 CRAN (R 3.3.2)
    ##  Rsamtools            * 1.26.1   2016-10-22 Bioconductor  
    ##  S4Vectors            * 0.12.2   2017-03-18 Bioconductor  
    ##  ShortRead            * 1.32.1   2017-03-18 Bioconductor  
    ##  stringi                1.1.3    2017-03-21 CRAN (R 3.3.2)
    ##  stringr                1.2.0    2017-02-18 CRAN (R 3.3.2)
    ##  SummarizedExperiment * 1.4.0    2016-10-18 Bioconductor  
    ##  withr                  1.0.2    2016-06-20 CRAN (R 3.3.0)
    ##  XVector              * 0.14.1   2017-03-18 Bioconductor  
    ##  yaml                   2.1.14   2016-11-12 CRAN (R 3.3.2)
    ##  zlibbioc               1.20.0   2016-10-18 Bioconductor
