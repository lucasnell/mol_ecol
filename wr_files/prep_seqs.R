
if(!'.wr_env' %in% ls(all.names = TRUE)) source('./wr_files/preamble.R')



sourceCpp('./wr_files/prep_seqs.cpp', env = .wr_env)


prep_seqs <- function(dna_in, read_len = 100, barcodes = c('CTCC', 'TGCA')) {
    bc_len <- nchar(barcodes[1])
    # Length to cut fragment to is simply the read length minus the barcode length and
    # the restriction enzyme's overhang length
    cut_len <- as.integer(read_len - bc_len)
    # applying the C++ inner function to each string in the character vector
    prep_seqs <- .wr_env$all_rm(as.character(dna_in), cut_len, barcodes[1], barcodes[2])
    # Adding names
    names(prep_seqs) <- paste0('seq_', rep(1:length(dna_in), 2), 
                               rep(c('_f', '_r'), each = length(dna_in)))
    dna_out <- .dna(prep_seqs)
    return(dna_out)
}

