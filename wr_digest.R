
# These are objects necessary for digesting a genome using restriction enzymes.

# If wr_preamble.R file hasn't been sourced, it needs to be.
# wr_preamble.R creates the .preamble_sourced object when it's run
if(!'.preamble_sourced' %in% ls(all.names = TRUE)) source('wr_preamble.R')


# Creating environment for these functions
digest_env <- new.env()

# This object stores the restriction enzymes I chose to test:
digest_env$chosen_enz <- c('ApeKI', 'BstBI', 'BspEI')

sourceCpp(code = 
'#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

vector<string> cpp_cut(string s, string prime5, string prime3) {
    
    string delimiter = prime5 + prime3;
    size_t last = 0;
    size_t next = 0;
    string token;
    vector<string> out_vec;
    // The first iteration needs to be treated differently...
    next = s.find(delimiter, last);
    if (next == string::npos) {
        out_vec.push_back(s);
        return out_vec;
    } else {
        token = s.substr(last, next - last + prime5.length());
        out_vec.push_back(token);
        last = next + prime5.length();
    }
    // Now I move onto more iterations...
    while ((next = s.find(delimiter, last + prime3.length())) != string::npos) {
        token = s.substr(last, next - last + prime5.length());
        out_vec.push_back(token);
        last = next + prime5.length();
    }
    out_vec.push_back(s.substr(last));
    
    return out_vec;
}

vector<string> cpp_cut_vec(vector<string> s_v, string prime5, string prime3) {
    
    vector<string> out_vec;
    vector<string> tmp_vec;
    string svi;
    for (int i = 0; i < s_v.size(); i++) {
        svi = s_v[i];
        tmp_vec = cpp_cut(svi, prime5, prime3);
        out_vec.insert(out_vec.end(), tmp_vec.begin(), tmp_vec.end());
    }
    
    return out_vec;
}

// [[Rcpp::export]]
vector<string> cpp_digest(CharacterVector DNAseq, CharacterVector cut_sites) {
    
    vector<string> out_vec;
    vector<string> tmp_vec;
    int cut_site_len = cut_sites.length();
    string dnasi, cs3, cs5;
    if ((cut_site_len % 2) != 0) {
        stop("cut_sites should have an even length");
    }
    for (int i = 0; i < DNAseq.length(); i++) {
        dnasi = DNAseq[i];
        cs5 = cut_sites[0];
        cs3 = cut_sites[1];
        tmp_vec = cpp_cut(dnasi, cs5, cs3);
        for (int j = 2; j < cut_site_len; j += 2) {
            cs5 = cut_sites[j];
            cs3 = cut_sites[j+1];
            tmp_vec = cpp_cut_vec(tmp_vec, cs5, cs3);
        }
        out_vec.insert(out_vec.end(), tmp_vec.begin(), tmp_vec.end());
    }
    return out_vec;
}
', env = digest_env)



# This function digests a genome (as a DNAStringSet object) using one of the 
# restriction enzymed used previously
digest_genome <- function(dna_ss, enzyme_name) {
    
    enz_df <- data_frame(
        enzyme = c('ApeKI', 'SbfI', 'PstI', 'EcoT22I', 'BstBI', 'AscI', 
                   'BspEI', 'AclI', 'FspI', 'MluI-HF', 'NruI-HF'),
        sites = list(c('G', 'CAGC', 'G', 'CTGC'), c('CCTGCA', 'GG'),
                     c('CTGCA', 'G'), c('ATGCA', 'T'), c('TT', 'CGAA'), 
                     c('GG', 'CGCGCC'), c('T', 'CCGGA'), c('AA', 'CGTT'),
                     c('TGC', 'GCA'), c('A', 'CGCGT'), c('TCG','CGA'))
    )
    
    enzyme_sites <- enz_df$sites[enz_df$enzyme == enzyme_name][[1]]
    
    # Converting from DNAStringSet to character vector and running C++ digestion
    dig_char <- digest_env$cpp_digest(as.character(dna_ss), enzyme_sites)
    
    dig_ss <- DNAStringSet(dig_char)
    
    return(dig_ss)
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
    } else if (class(dna_ss) != 'list' | any(lapply(dna_ss, class) != 'DNAStringSet')) {
        stop('dna_ss must be a list or DNAStringSet object')
    }
    if (is.null(file_names)) {
        file_names <- paste0('frags_', 1:length(dna_ss), '.fa.gz')
    } else if (length(file_names) != length(dna_ss)) {
        stop('dna_ss and file_names must be same length')
    }
    # Using purrr::map and magrittr::set_names to set sequence names before writing
    # (Otherwise, each sequence-name line is just ">")
    write_list <- map(dna_ss, ~ magrittr::set_names(.x, paste0('seq_', 1:length(.x))))
    for (i in 1:length(write_list)) {
        writeFasta(write_list[[i]], file = file_names[i], mode = 'w', 
                   compress = grepl('.gz', file_names[i]))
        cat(file_names[i], 'written... \n')
    }
    return(invisible(NULL))
}



# # ---------------
# # Example usage:
# # ---------------
# 
# set.seed(9)
# test_ <- read_fasta('./genome_data/aphid_genome.fa.gz')
# test_ <- test_[sample(1:length(test_), 1e4L)]
# enzyme <- 'ApeKI'
# 
# digest <- digest_genome(test_, 'ApeKI')
# 
# 
# write_fastas(digest, '~/Desktop/digestion.fa.gz')

