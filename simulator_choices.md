Choosing the sequencing simulator
================
Lucas Nell
19 March 2017

-   [Reading Excel sheets](#reading-excel-sheets)

*Updated 24 March 2017*

Tables are from the following paper (link [here](http://www.nature.com/nrg/journal/v17/n8/full/nrg.2016.57.html)):

Escalona, M., S. Rocha, and D. Posada. 2016. A comparison of tools for the simulation of genomic next-generation sequencing data. *Nature Reviews Genetics* **17**:459–469.

**Loading packages:**

``` r
suppressPackageStartupMessages({
    library(dplyr)
    library(readxl)
})
```

Reading Excel sheets
--------------------

``` r
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
```

    ##      abbrev                                                 full
    ## 1       454                                   454 pyrosequencing
    ## 2        AL                                           alignments
    ## 3       API                    application programming interface
    ## 4        AU                                    academic use only
    ## 5       BSD                       Berkeley software distribution
    ## 6    CCANCL  creative commons attribution non-commercial license
    ## 7       CCS                        circular consensus sequencing
    ## 8       CLI                               command line interface
    ## 9       CLR                                 continuous long read
    ## 10     CNVs                                 copy number variants
    ## 11       DF                                      default profile
    ## 12       FA                                                FASTA
    ## 13       FO                                               format
    ## 14       FQ                                                FASTQ
    ## 15        G                                             genomics
    ## 16      GPL                               general public license
    ## 17       GU                           guide to generate profiles
    ## 18      GUI                             graphical user interface
    ## 19       GV                                     genomic variants
    ## 20   indels                          insertions and/or deletions
    ## 21     INVs                                           inversions
    ## 22     LGPL                        lesser general public license
    ## 23        M                                         metagenomics
    ## 24      MGC                                metagenomic community
    ## 25      MIT                Massachusetts Institute of Technology
    ## 26       MP                                            mate pair
    ## 27       NA                                       not applicable
    ## 28       NA                                       not applicable
    ## 29 Nanopore                         Oxford Nanopore Technologies
    ## 30      NGS                          next- generation sequencing
    ## 31      NGS                           next-generation sequencing
    ## 32       NP                               no parallel processing
    ## 33        P     parallel processing that accepts multi-threading
    ## 34       PA                                            parameter
    ## 35   PacBio                                  Pacific Biosciences
    ## 36       PE                                           paired end
    ## 37      PLO                                               ploidy
    ## 38       PR                                              profile
    ## 39      PRO                                 proprietary software
    ## 40       QS                                        quality score
    ## 41       RE                                                reads
    ## 42  Ref seq                                   reference sequence
    ## 43   Sanger                                    Sanger sequencing
    ## 44       SE                                           single end
    ## 45      SFF                             standard flowgram format
    ## 46     SNPs                      single-nucleotide polymorphisms
    ## 47    SOLiD sequencing by oligonucleotide ligation and detection
    ## 48     STRs                                 short tandem repeats
    ## 49       SW                specific software to generate profile
    ## 50      TRA                                        translocation

``` r
all_tables %>% 
    filter(
        grepl('Illumina', Technology),
        grepl('SE', Run_types),
        PCR == 'Yes',
        grepl('Mac', Operating_system),
        Open_source == 'Yes',
        SNPs == 'Yes'
        ) %>% 
    # select(-Technology, -Run_types, -PCR, -Operating_system, -Open_source, -SNPs) %>%
    # select(Simulators) %>%
    as.data.frame
```

    ##   Simulators               Technology  G_vs_M     Run_types PCR  QS out_RE
    ## 1    Grinder 454, Illumina and Sanger G and M SE, PE and MP Yes Yes    Yes
    ##   AL FO Programming_language          Operating_system Processing
    ## 1 No FQ                 Perl Windows, Linux and Mac OS          P
    ##   Open_source SNPs
    ## 1         Yes  Yes

``` r
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
```

    ##   Simulators               Technology  G_vs_M     Run_types PCR  QS out_RE
    ## 1        ART  454, Illumina and SOLiD       G SE, PE and MP Yes Yes    Yes
    ## 2    Flowsim                      454       G     SE and PE Yes Yes    Yes
    ## 3    Grinder 454, Illumina and Sanger G and M SE, PE and MP Yes Yes    Yes
    ##    AL         FO Programming_language          Operating_system Processing
    ## 1 Yes SFF and FQ         C++ and Perl Windows, Linux and Mac OS          P
    ## 2  No        SFF              Haskell                     Linux          P
    ## 3  No         FQ                 Perl Windows, Linux and Mac OS          P
    ##   Open_source SNPs
    ## 1         Yes <NA>
    ## 2         Yes <NA>
    ## 3         Yes  Yes

``` r
with(all_tables, filter(all_tables, eval(eval_list[['snps']]))) %>% as.data.frame
```

    ##               Simulators                             Technology  G_vs_M
    ## 1  DWGSIM (DNA analysis)         Illumina, SOLiD and IonTorrent       G
    ## 2                  EAGLE   454, Illumina, PacBio and IonTorrent       G
    ## 3               FASTQSim Illumina, SOLiD, PacBio and IonTorrent G and M
    ## 4                 GemSim                       454 and Illumina G and M
    ## 5                Grinder               454, Illumina and Sanger G and M
    ## 6                  Mason               454, Illumina and Sanger       G
    ## 7                   pIRS                               Illumina G and M
    ## 8                ReadSim                    PacBio and Nanopore       G
    ## 9                   SInC                               Illumina       G
    ## 10                 wgsim                     Illumina and SOLiD       G
    ##        Run_types PCR  QS out_RE  AL        FO Programming_language
    ## 1  SE, PE and MP  No Yes    Yes  No        FQ   C, Perl and Python
    ## 2      SE and PE  No Yes    Yes Yes        FQ                  C++
    ## 3             SE  No Yes    Yes  No        FQ      Bash and Python
    ## 4      SE and PE  No Yes    Yes  No        FQ               Python
    ## 5  SE, PE and MP Yes Yes    Yes  No        FQ                 Perl
    ## 6  SE, PE and MP  No Yes    Yes Yes FA and FQ                  C++
    ## 7             PE  No Yes    Yes  No        FQ         C++ and Perl
    ## 8             SE  No Yes    Yes  No        FQ               Python
    ## 9             PE  No Yes    Yes  No        FQ                  C++
    ## 10            SE  No Yes    Yes  No        FQ                    C
    ##             Operating_system Processing Open_source SNPs
    ## 1                      Linux          P         Yes  Yes
    ## 2                      Linux   NP and P         Yes  Yes
    ## 3                      Linux   NP and P         Yes  Yes
    ## 4  Windows, Linux and Mac OS          P         Yes  Yes
    ## 5  Windows, Linux and Mac OS          P         Yes  Yes
    ## 6  Windows, Linux and Mac OS          P         Yes  Yes
    ## 7                      Linux   NP and P         Yes  Yes
    ## 8  Windows, Linux and Mac OS          P         Yes  Yes
    ## 9                      Linux   NP and P          No  Yes
    ## 10                     Linux          P         Yes  Yes
