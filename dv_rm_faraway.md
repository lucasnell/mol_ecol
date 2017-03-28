Remove areas faraway from cut sites
================
Lucas Nell
28 March 2017

*Updated 28 March 2017*

This script provides background for the function that removes regions from digested fragments that will not be sequenced, those that are far from restriction enzyme cut sites

**Loading packages:**

``` r
suppressPackageStartupMessages({
    library(purrr)
    library(dplyr)
    library(ShortRead)
})
```

Davey, J. W., P. A. Hohenlohe, P. D. Etter, J. Q. Boone, J. M. Catchen, and M. L. Blaxter. 2011. Genome-wide genetic marker discovery and genotyping using next-generation sequencing. *Nature Reviews Genetics* **12**:499–510. Available from <http://dx.doi.org/10.1038/nrg3012>

Elshire, R. J., J. C. Glaubitz, Q. Sun, J. A. Poland, K. Kawamoto, E. S. Buckler, and S. E. Mitchell. 2011. A robust, simple genotyping-by-sequencing (GBS) approach for high diversity species. *PLOS ONE* **6**:e19379. Available from <http://dx.plos.org/10.1371/journal.pone.0019379>
