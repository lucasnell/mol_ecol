
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


// [[Rcpp::export]]
vector<string> cpp_cut(string s, string prime5, string prime3) {
    
    string delimiter = prime5 + prime3;
    size_t last = 0;
    size_t next = 0;
    string token;
    vector<string> out_vec;
    next = s.find(delimiter, last);
    if (next == string::npos) {
        out_vec.push_back(s);
        return out_vec;
    } else {
        token = s.substr(last, next - last + prime5.length());
        out_vec.push_back(token);
        last = next + prime5.length();
    }
    while ((next = s.find(delimiter, last + prime3.length())) != string::npos) {
        token = s.substr(last, next - last + prime5.length());
        out_vec.push_back(token);
        last = next + prime5.length();
    }
    out_vec.push_back(s.substr(last));
    
    return out_vec;
}

// [[Rcpp::export]]
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





string cpp_merge_str(CharacterVector in_strings) {
    
    int num_char = in_strings.size();
    string out_str;
    
    for(int j=0; j < num_char; j++) {
        out_str += in_strings[j];
    }
    
    return out_str;
}


// [[Rcpp::export]]
vector<string> cpp_is_digest(CharacterVector DNAseq, CharacterVector cut_sites) {
    
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



test_fun <- function(input_vec, delim){
    .z <- strsplit(input_vec, paste(delim, collapse = ''))
    for (i in 1:length(.z)){
        .x <- .z[[i]]
        if (length(.x) == 1) {
            next
        }
        .x[1] <- paste0(.x[1], delim[1])
        if (length(.x) > 1) {
            .x[length(.x)] <- paste0(delim[2], .x[length(.x)])
            if (length(.x) > 2) {
                .x[2:(length(.x)-1)] <- paste0(delim[2], .x[2:(length(.x)-1)], delim[1])
            }
        }
        .z[[i]] <- .x
    }
    return(unlist(.z))
}


seqs <- as.character(test_[1:100])
delim <- c('G', 'CAGC')

test_fun(seqs[2], delim) %>% head

x <- test_fun(seqs, delim)
y <- digest_env$cpp_cut_vec(seqs, delim[1], delim[2])
identical(x, y)
str(x); str(y)
# which(x != y)

for (i in 1:length(x)) {
    if (!identical(x[i], y[i])){
        print(i)
        break
    }
}

x[3]; y[3]

for (i in 61:63) {
    print(x[i]); print(y[i])
    cat('\n')
}; rm(i)


digest_env$cpp_cut_vec(c('ABCD', 'AABCDEF'), 'B', 'C')

digest_env$cpp_cut_vec("GCGTGGCAGCAGCAAATTAGCCACA", 'G', 'CAGC')
test_fun("GCGTGGCAGCAGCAAATTAGCCACA", c('G', 'CAGC'))





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
    dna_char <- as.character(dna_ss)
    
    for (i in seq(1, length(enzyme_sites), 2)) {
        call_list <- as.list(c(DNAseq = dna_char, verbose = FALSE, enzyme_sites))
        dig <- do.call(insilico.digest, call_list)
    }
    
    dig_ss <- DNAStringSet(dig)
    
    return(dig_ss)
}


insilico.digest_ <- function (DNAseq, recognition_code, cut_site_5prime, 
                               cut_site_3prime) {
    frag1 <- lapply(DNAseq, digest_env$cpp_cut, prime5 = cut_site_5prime, 
                    prime3 = cut_site_3prime)
    frag2 <- c(frag1, recursive = TRUE)
    return(frag2)
}

r_is <- function(DNAseq, cut_sites) {
    rc <- paste0(cut_sites[1:2], collapse = '')
    frag1 <- insilico.digest_(DNAseq, rc, cut_sites[1], cut_sites[2])
    if (length(cut_sites) == 4) {
        rc <- paste0(cut_sites[3:4], collapse = '')
        frag2 <- insilico.digest_(frag1, rc, cut_sites[3], cut_sites[4])
        return(frag2)
    }
    return(frag1)
}


set.seed(9)
test_ <- read_fasta('./genome_data/aphid_genome.fa.gz')
test_ <- test_[sample(1:length(test_), 1e3L)]
cut_sites <- c('G', 'CAGC', 'G', 'CTGC')


# system.time(dig1 <- SimRAD::insilico.digest(as.character(test_), cut_sites[1], cut_sites[2], verbose = FALSE))
system.time(dig1 <- SimRAD::insilico.digest.internal(as.character(test_), paste0(cut_sites[1], cut_sites[2]), cut_sites[1], cut_sites[2]))
system.time(dig2 <- digest_env$cpp_is_digest(as.character(test_), cut_sites[1:2]))
dig2 <- dig2[dig2 != 'CAGC']
# system.time(dig3 <- r_is(as.character(test_), cut_sites[1:2]))
identical(sort(dig1), sort(dig2))
identical(nchar(dig1), nchar(dig2))

length(dig1)
length(dig2)
# length(dig3)


cuts <- stringr::str_count(as.character(test_), paste0(cut_sites[1], cut_sites[2]))
sum(cuts) + length(cuts)

for (i in 1:length(dig1)) {
    if (!identical(dig1[i], dig2[i]) & i != 1250) {
        print(i)
        break
    }
}
for (i in 1249:1251) {
    print(dig1[i]); print(dig2[i])
    cat('\n\n')
}

dig1[1250] %>% substr(nchar(.) - 10, nchar(.))
dig2[1250] %>% substr(nchar(.) - 10, nchar(.))

dig1[1250] %>% nchar
dig2[1250] %>% nchar

dig1[1251] %>% substr(1, 10)
dig2[1251] %>% substr(1, 10)

which(grepl("CAGCAGTCAGACGGATGTTCGCAGTAACCGCGCTATAACCAGAAATGCTGACATTGATTGGCCACCCCTATCCCCACACCAATCCGTCTGTGATTTTTTTTATGGTGACATCTGAAAAATGTCGTGTGCCACACTCGCCCCACAACCATTTCCGATCTCGACTTGGGGTTCGCCGAAAAAATCAATAATATTCCACCTGATACCTTGCACCAGGCATTGAACATTTTCAAAAATCGACTTGACAAATGCATTCGCAGAAATGGGCCGATGACTGATATCGTTTAAAAAAAAAAAAATATTATTTTAAAATATTGTAGAATACATCGATTTTAAAAAAATGCTTTCGTCCTTGGTATTACATTTAAAACTGTTTTTACATTGACACGGCAAAAAGTATTTAATGTTATAACTTATATTATATGAGCATGTTACATGTGTATACCGTGTACGTATAGTGCAAAACCTATTACAGGCTACAGTTAAAATACACAAAAAATTACCGAACACAATATGTCGAATAAAGTATGATTTATCGGGGGATGGAGTGGTCGAGCGGACTAAGGCGTCGGTTGCGACGCACACCGACGCCGGTTCAAAGCTTCGGCCGGGGGCGGCATTTTTCTTCGGGCAAGTCACGGTAGCCGGAGAGAACCCCTACTCGGGCATGGAAGATACCTATGGGTGCCCAATAAAAAATCTGCCAAACCACACGTGTTTTACCACTTCCCTACAAAAACAAAAAACAGCTAATGGCCATAGTTGCCGGGCTTTACGATCAATTAAAAAAAAAAAAAAGTATGATTTGTTATTATGATTTATAATAATATACATACCTCAACTAAAAATTCAAGAGACAATAGGTATGATAATATCTTGGTAAAGGGTATTGGTACCTATTGGGTATATCGTATAATATAAAATGTATTGGTTGTAATATAATTCCGTAAAATACGTATTTTAATAATTTAATTTATGTTCCTGTAACCTGCGAAATT", as.character(test_)))


stringr::str_locate(test_[220])

# rm(dig1, dig2, test_, cs5, cs3, rc)



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
