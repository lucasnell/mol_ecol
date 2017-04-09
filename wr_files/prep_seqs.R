
if(!'.wr_env' %in% ls(all.names = TRUE)) source('./wr_files/preamble.R')



sourceCpp('./wr_files/prep_seqs.cpp', env = .wr_env)



# Preps one set of sequences for simulation
.wr_env$prep_one_seqs <- function(dna_in, samp_num, f_len, r_len, paired, barcodes) {
    if (paired) {
        prep_seqs <- .wr_env$prep_p(as.character(dna_in), f_len, r_len,
                                    barcodes[1], barcodes[2])
        names(prep_seqs) <- paste0('samp_', samp_num, '_seq_', 1:length(dna_in))
    } else {
        prep_seqs <- .wr_env$prep_s(as.character(dna_in), f_len, r_len, 
                                    barcodes[1], barcodes[2])
        names(prep_seqs) <- paste0('samp_', samp_num, '_seq_', rep(1:length(dna_in), 2), 
                                   rep(c('_f', '_r'), each = length(dna_in)))
    }
    dna_out <- .dna(prep_seqs)
    return(dna_out)
}

prep_seqs <- function(dna_list_in, read_len = 100, paired = TRUE, 
                      barcodes = c('CTCC', 'TGCA')) {
    f_len <- read_len - nchar(barcodes[1])
    r_len <- read_len - nchar(barcodes[2])
    dna_list_tmp <- lapply(1:length(dna_list_in), 
                           function(.i) {
                               .wr_env$prep_one_seqs(dna_list_in[[.i]], .i,
                                                     f_len, r_len, paired, barcodes)
                           })
    dna_list_out <- .dna_list(dna_list_tmp)
    return(dna_list_out)
}

