
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
    # Adding names
    names(cut_dna_ss) <- paste0('seq_', rep(1:length(dna_ss), 2), 
                                rep(c('_f', '_r'), each = length(dna_ss)))
    return(cut_dna_ss)
}
