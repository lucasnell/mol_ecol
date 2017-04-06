sourceCpp(code =
'#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// Function to reverse a string
string reverse_str(string in_str) {
    string out_str = in_str;
    int n = in_str.length();
    // Swap character starting from two corners
    for (int i=0; i<n/2; i++) {
        swap(out_str[i], out_str[n-i-1]);
    }
    return out_str;
}

// Function to return complement of a string
// [[Rcpp::export]]
string comp_str(string in_str) {
    string out_str = in_str;
    // CharacterVector nts = CharacterVector::create("A", "C", "G", "T");
    // CharacterVector repl = CharacterVector::create("T", "G", "C", "A");
    // // Switch characters
    // for (int i=0; i < 4; i++) {
    //     replace(out_str.begin(), out_str.end(), nts[i], repl[i]);
    // }
    replace(out_str.begin(), out_str.end(), "C", "G");
    return out_str;
}

// // [[Rcpp::export]]
// vector<string> rcomp(vector<string> in_strings) {
// 
//     int num_char = in_strings.size();
//     string rev_str, rcomp_str;
//     vector<string> out_strings(num_char);
//     
//     for(int i=0; i < num_char; i++) {
//         rev_str = reverse_str(in_strings[i]);
//         rcomp_str = comp_str(rev_str);
//         out_strings[i] = rcomp_str;
//     }
//     
//     return out_strings;
// }
')

cppFunction(
'
std::string cpp_merge_str(CharacterVector in_strings) {
    int num_char = in_strings.size();
    std::string out_str;
    for(int j=0; j < num_char; j++) {
        out_str += in_strings[j];
    }
    return out_str;
}
')


cppFunction(
"std::string comp_str(std::string in_str) {
    std::string out_str = in_str;
    // char nts[] =  {'A', 'C', 'G', 'T'};
    // char repl[] = {'T', 'G', 'C', 'A'};
    // // Switch characters
    // for (int i=0; i < 4; i++) {
    //     std::replace(out_str.begin(), out_str.end(), nts[i], repl[i]);
    // }
    std::replace(out_str.begin(), out_str.end(), 'A', 'T');
    std::replace(out_str.begin(), out_str.end(), 'C', 'G');
    // std::replace(out_str.begin(), out_str.end(), 'G', 'C');
    // std::replace(out_str.begin(), out_str.end(), 'T', 'A');
    // std::replace(out_str.begin(), out_str.end(), 'A', 'T');
    return out_str;
}
")


comp_str('AAAACCCGGGTTT')


# Left off --> It's changing all As that get turned to Ts back into As at the end...


# Inner function that removes areas far from cut sites from one sequence
cppFunction(
'CharacterMatrix all_rm(CharacterVector all_s, int cut_len) {
    std::string f_strand, r_strand, s;
    int seq_len, vec_len = all_s.length(); // = s.size();
    CharacterMatrix out_mat(vec_len,2);
    for (int i = 0; i < vec_len; i++) {
        s = all_s[i];
        seq_len = s.size();
        if (seq_len < cut_len) {
            f_strand = s;
            r_strand = s;
        } else {
            f_strand = s.substr(0, cut_len);
            r_strand = s.substr(seq_len - cut_len, cut_len);
        }
        out_mat(i,0) = f_strand;
        out_mat(i,1) = r_strand;
    }
    return out_mat;
}
')


prep_seqs <- function(dna_ss, read_len = 100, barcodes = c('CTCC', 'TGCA')) {
    bc_len <- nchar(barcodes[1])
    # Length to cut fragment to is simply the read length minus the barcode length and
    # the restriction enzyme's overhang length
    cut_len <- as.integer(read_len - bc_len)
    # Extracting character vector from input DNAStringSet
    character_seqs <- as.character(dna_ss)
    # applying that inner function to each string in the character vector
    cut_seq_char <- all_rm(character_seqs, cut_len)
    # Creating DNAStringSet objects for forward and reverse strands.
    # Because I am doing unidirectional (forward-only) sequencing, I need to make the 
    # reverse strands reverse complements of the forward one.
    f_strands <- DNAStringSet(paste0(barcodes[1], cut_seq_char[,1]))
    r_strands <- reverseComplement(DNAStringSet(paste0(barcodes[2], cut_seq_char[,2])))
    # Now I combine them
    cut_dna_ss <- append(f_strands, r_strands)
    # Adding names
    names(cut_dna_ss) <- paste0('seq_', rep(1:length(dna_ss), 2), 
                                rep(c('_f', '_r'), each = length(dna_ss)))
    return(cut_dna_ss)
}






prep_seqs2 <- function(dna_ss, read_len = 100, barcodes = c('CTCC', 'TGCA')) {
    bc_len <- nchar(barcodes[1])
    # Length to cut fragment to is simply the read length minus the barcode length and
    # the restriction enzyme's overhang length
    cut_len <- as.integer(read_len - bc_len)
    # Extracting character vector from input DNAStringSet
    character_seqs <- as.character(dna_ss)
    # applying that inner function to each string in the character vector
    cut_seq_char <- all_rm(character_seqs, cut_len)
    # Creating DNAStringSet objects for forward and reverse strands.
    # Because I am doing unidirectional (forward-only) sequencing, I need to make the 
    # reverse strands reverse complements of the forward one.
    f_strands <- DNAStringSet(paste0(barcodes[1], cut_seq_char[,1]))
    r_strands <- reverseComplement(DNAStringSet(paste0(barcodes[2], cut_seq_char[,2])))
    # Now I combine them
    cut_dna_ss <- append(f_strands, r_strands)
    # Adding names
    names(cut_dna_ss) <- paste0('seq_', rep(1:length(dna_ss), 2), 
                                rep(c('_f', '_r'), each = length(dna_ss)))
    return(cut_dna_ss)
}
