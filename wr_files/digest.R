
# These are objects necessary for digesting a genome using restriction enzymes.

# If wr_preamble.R file hasn't been sourced, it needs to be.
# wr_preamble.R creates the .preamble_sourced object when it's run
if(!'.wr_env' %in% ls(all.names = TRUE)) source('./wr_files/preamble.R')


# This object stores the restriction enzymes I chose to test:
.wr_env$chosen_enz <- c('ApeKI', 'BstBI', 'BspEI')

sourceCpp('./wr_files/digest.cpp', env = .wr_env)



# This function digests a genome (as a dna object) using one of the 
# restriction enzymed used previously
digest_genome <- function(dna, enzyme_name) {
    
    enz_df <- data_frame(
        enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
                   'BspEI', 'AclI', 'FspI', 'MluI-HF', 'NruI-HF'),
        sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                     c('CTGCA', 'G'), c('ATGCA', 'T'), c('TT', 'CGAA'), 
                     c('GG', 'CGCGCC'), c('T', 'CCGGA'), c('AA', 'CGTT'),
                     c('TGC', 'GCA'), c('A', 'CGCGT'), c('TCG','CGA'))
    )
    
    enzyme_sites <- enz_df$sites[enz_df$enzyme == enzyme_name][[1]]
    
    # Running C++ digestion
    dig_char <- .wr_env$cpp_digest(dna, enzyme_sites)
    
    dig_dna <- .dna(dig_char)
    
    return(dig_dna)
}



# # ---------------
# # Example usage:
# # ---------------
# 
# set.seed(9)
# test_ <- read_fasta('./genome_data/aphid_genome.fa.gz')
# test_ <- test_[sample(1:length(test_), 1e2L)]
# enzyme <- 'ApeKI'
# 
# digest <- digest_genome(test_, 'ApeKI')
# 
# write_fastas(digest, '~/Desktop/digestion.fa.gz')

