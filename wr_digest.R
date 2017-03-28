
# These are objects necessary for digesting a genome using restriction enzymes.


# This function digests a genome (as a DNAStringSet object) using one of the 
# restriction enzymed used previously
digest_genome <- function(enzyme_name, dna_ss) {
    
    enz_df <- data_frame(
        enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
                   'BspEI', 'AclI', 'FspI', 'MluI-HF', 'NruI-HF'),
        sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                     c('CTGCA', 'G'), c('ATGCA', 'T'), c('TT', 'CGAA'), 
                     c('GG', 'CGCGCC'), c('T', 'CCGGA'), c('AA', 'CGTT'),
                     c('TGC', 'GCA'), c('A', 'CGCGT'), c('TCG','CGA'))
    )
    
    enzyme_sites <- enz_df$sites[enz_df$enzyme == enzyme_name][[1]]
    
    names(enzyme_sites) <- 
        c('cut_site_5prime1', 'cut_site_3prime1', 
          'cut_site_5prime2', 'cut_site_3prime2',  
          'cut_site_5prime3', 'cut_site_3prime3', 
          'cut_site_5prime4', 'cut_site_3prime4')[1:length(enzyme_sites)]
    
    # Converting from DNAStringSet to character vector
    dna_str <- paste(dna_ss, collapse = '')
    
    call_list <- as.list(c(DNAseq = dna_str, verbose = FALSE, enzyme_sites))
    
    dig <- do.call(insilico.digest, call_list)
    
    return(dig)
}


# This reads a fasta file into a DNAStringSet object
read_fasta <- function(file_name) {
    fasta <- readFasta(file_name)
    dna_ss <- sread(fasta)
    return(dna_ss)
}



# This function can write ≥ 1 DNAStringSet objects wrapped in a list to a fasta file.
# You can also input a single DNAStringSet
write_fastas <- function(dna_ss, file_names = NULL) {
    if (class(dna_ss) == 'DNAStringSet') {
        dna_ss <- list(dna_ss)
    } else if (class(dna_ss) != 'list') {
        stop('dna_ss must be a list or DNAStringSet object')
    }
    if (is.null(file_names)) {
        file_names <- paste0('frags_', 1:length(dna_ss), '.fa.gz')
    } else if (length(file_names) != length(dna_ss)) {
        stop('dna_ss and file_names must be same length')
    }
    # Using map and set_names from purrr to set sequence names before writing
    # (Otherwise, each sequence-name line is just ">")
    write_list <- map(dna_ss, ~ purrr::set_names(.x, paste0('seq_', 1:length(.x))))
    for (i in 1:length(write_list)) {
        writeFasta(write_list[[i]], file = file_names[i], mode = 'w', 
                   compress = grepl('.gz', file_names[i]))
        cat(file_names[i], 'written... \n')
    }
}
