#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]






using namespace Rcpp;
using namespace std;


/*
 
 #include <RcppParallel.h>
 
 // [[Rcpp::depends(RcppParallel)]]
 
 using namespace RcppParallel;
 
 
 */



/*
 ------------
 Shuffle integer vector
 ------------
*/

IntegerVector shuffle_int(IntegerVector in_vec) {
    random_shuffle(in_vec.begin(), in_vec.end());
    return in_vec;
}




/*
 ------------
 Make sequence from a vector of frequencies
 ------------
*/

// [[Rcpp::export]]
CharacterVector cpp_make_seq(IntegerVector acgt) {
    int total_nt = sum(acgt);
    CharacterVector bases = CharacterVector::create('A', 'C', 'G', 'T');
    string nt_i;
    CharacterVector seq_out(total_nt);
    int k=0;
    
    for (int i = 0; i < 4; i++) {
        nt_i = bases[i];
        for (int j = 0; j < acgt[i]; j++) {
            seq_out[k] = nt_i;
            k+=1;
        }
    }
    
    return seq_out;
}





/*
 ------------
 Merge a vector of strings into one
------------
*/

// [[Rcpp::export]]
string cpp_merge_str(CharacterVector in_strings) {
    
    int num_char = in_strings.size();
    string out_str;
    
    for(int j=0; j < num_char; j++) {
        out_str += in_strings[j];
    }
    
    return out_str;
}





/*
 ------------
 Split one string into a vector
------------
*/

CharacterVector cpp_str_split1(string in_string, int n = 1) {
    
    int num_substr = in_string.length() / n;
    
    CharacterVector out(num_substr);
    
    for(int j=0; j < num_substr; j++) {
        out[j] = in_string.substr(j*n, n);
    }
    
    return out;
}



/*
 ------------
 RcppArmadillo's version of R's sample function
------------
*/

IntegerVector cpp_sample(IntegerVector x, int size,
                         bool replace = false, 
                         NumericVector prob = NumericVector::create()) {
    IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
    return ret;
}


/*
 ------------
 Random integer generation, within a range
 ------------
*/

IntegerVector int_sampler(int num_samps, float range_min, float range_max) {
    IntegerVector out_vec;
    NumericVector tmp_vec = runif(num_samps, range_min - 1, range_max);
    out_vec = ceil(tmp_vec);
    return out_vec;
}


/*
 ------------
 Random integer generation, within a range, NO replacement
------------
*/

IntegerVector int_sampler_nr(int num_samps, float range_min, float range_max) {
    if ((range_max - range_min + 1) < (num_samps)) {
        stop("n_samps > range of values; can't sample without replacement");
    }
    IntegerVector out_vec;
    IntegerVector uniques;
    NumericVector tmp_vec = runif(num_samps, range_min - 1, range_max);
    out_vec = ceil(tmp_vec);
    uniques = unique(out_vec);
    int k = 0;
    while (out_vec.size() != uniques.size()) {
        tmp_vec = runif(num_samps, range_min - 1, range_max);
        out_vec = ceil(tmp_vec);
        uniques = unique(out_vec);
        k += 1;
        if (k > 10) {
            out_vec = cpp_sample(Range(range_min, range_max), num_samps);
            break;
        }
    }
    return out_vec;
}

/*
 ------------
 Get the site locations for 1 sequence
 ------------
*/

// [[Rcpp::export]]
IntegerVector cpp_one_sites(IntegerVector freq_len_row) {
    
    int samp_n = freq_len_row[0];
    if (samp_n == 0) {
        return IntegerVector::create();
    }
    int seq_length = freq_len_row[1];
    IntegerVector out_locs = int_sampler_nr(samp_n, 1, seq_length);
    return out_locs;
}





/*
 ------------
 Get the site locations for multiple sequences
 ------------
*/

// [[Rcpp::export]]
IntegerVector cpp_get_sites(IntegerMatrix freq_len) {
    
    int n = freq_len.nrow();
    int total_seqs = sum(freq_len(_,0));
    int samp_n;
    int j = 0;
    
    IntegerVector ran_locs(total_seqs);
    IntegerVector tmp_locs;
    
    for(int i = 0; i < n; i++) {
        samp_n = freq_len(i,0);
        if (samp_n == 0) {
            continue;
        }
        tmp_locs = cpp_one_sites(freq_len(i,_));
        ran_locs[Range(j, j + samp_n - 1)] = tmp_locs;
        j += samp_n;
    }
    
    return ran_locs;
}









/*
 ------------
 Change sites for 1 sequence
 ------------
*/

// [[Rcpp::export]]
CharacterVector cpp_change(string seq, IntegerVector positions,
                           IntegerMatrix freq_mat, int n_samps) {
    
    int n_pos = positions.size();
    CharacterVector out_vec(n_samps);
    
    if (n_pos == 0) {
        for (int i = 0; i < n_samps; i++) {
            out_vec[i] = seq;
        }
        return out_vec;
    }
    
    CharacterMatrix seq_mat(n_samps,seq.size());
    CharacterVector seq_out(n_samps);
    
    CharacterVector seq_vec = cpp_str_split1(seq);
    for (int i = 0; i < n_samps; i++) {
        seq_mat(i,_) = seq_vec;
    }

    int ran_r = rand() % freq_mat.nrow();
    IntegerVector ran_c = Range(0, 4 - 1);
    random_shuffle(ran_c.begin(), ran_c.end());
    IntegerVector ran_row = freq_mat(ran_r, _);
    IntegerVector ran_freq = ran_row[ran_c];

    CharacterVector tmp_seq;

    for (int i = 0; i < n_pos; i++) {
        tmp_seq = cpp_make_seq(ran_freq);
        seq_mat(_,positions[i]) = tmp_seq;
    }
    for (int i = 0; i < n_samps; i++) {
        tmp_seq = seq_mat(i,_);
        out_vec[i] = cpp_merge_str(tmp_seq);
    }
    
    return out_vec;
}









// IntegerVector cpp_one_sites_par(RMatrix<int>::Row freq_len_row) {
//     
//     int samp_n = freq_len_row[0];
//     IntegerVector out_locs;
//     if (samp_n == 0) {
//         out_locs = IntegerVector::create();
//         return out_locs;
//     }
//     int seq_length = freq_len_row[1];
//     out_locs = int_sampler_nr(samp_n, 0, seq_length - 1);
//     return out_locs;
// }
// 
// CharacterVector cpp_changep(string seq, IntegerVector positions,
//                             RMatrix<int> freq_mat, int n_samps) {
//     
//     int n_pos = positions.size();
//     CharacterVector out_vec(n_samps);
//     CharacterVector seq_out(n_samps);
//     // RVector<string> out_vec(seq_out);
//     
//     if (n_pos == 0) {
//         for (int i = 0; i < n_samps; i++) {
//             out_vec[i] = seq;
//         }
//         return out_vec;
//     }
//     
//     CharacterMatrix seq_mat(n_samps,seq.size());
//     
//     CharacterVector seq_vec = cpp_str_split1(seq);
//     for (int i = 0; i < n_samps; i++) {
//         seq_mat(i,_) = seq_vec;
//     }
//     
//     int ran_r = rand() % freq_mat.nrow();
//     IntegerVector ran_c = Range(0, 4 - 1);
//     random_shuffle(ran_c.begin(), ran_c.end());
//     RMatrix<int>::Row ran_row = freq_mat.row(ran_r);
//     IntegerVector ran_freq;
//     for (int i=0; i < freq_mat.ncol(); i++) {
//         ran_freq[i] = ran_row[ran_c[i]];
//     }
//     
//     CharacterVector tmp_seq;
//     
//     for (int i = 0; i < n_pos; i++) {
//         tmp_seq = cpp_make_seq(ran_freq);
//         seq_mat(_,positions[i]) = tmp_seq;
//     }
//     for (int i = 0; i < n_samps; i++) {
//         tmp_seq = seq_mat(i,_);
//         out_vec[i] = cpp_merge_str(tmp_seq);
//     }
//     
//     return out_vec;
// }
// 
// /*
// ------------
// Change sites for multiple locations
// ------------
// */
// 
// 
// 
// struct new_seq : public Worker {
//     
//     // frequency-length matrix
//     const RMatrix<int> freq_len;
//     // sequence vector
//     const RVector<string> seqs;
//     // frequency matrix
//     const RMatrix<int> freq_mat;
//     // number of samples
//     const int n_samps;
//     
//     // destination vector
//     RVector<string> output;
//     
//     // initialize with source and destination
//     new_seq(const IntegerMatrix freq_len, const CharacterVector seqs, 
//             const IntegerMatrix freq_mat, const int n_samps,
//             CharacterVector output)
//         : freq_len(freq_len), seqs(seqs), freq_mat(freq_mat), n_samps(n_samps), 
//           output(output) {}
//     
//     // // Iteratively create new sequences
//     void operator()(size_t begin, size_t end) {
//         for (size_t i = begin; i < end; i++) {
//             RMatrix<int>::Row freq_len_row = freq_len.row(i);
//             string seq = seqs[i];
//             IntegerVector sites = cpp_one_sites_par(freq_len_row);
//             CharacterVector new_seqs = cpp_changep(seq, sites, freq_mat, n_samps);
//             for (int j = 0; j < new_seqs.size(); j++) {
//                 output[((i * 10) + j)] = new_seqs[j];
//             }
//         }
//     }
// };
// 
// 
// 
// 
// 
// // [[Rcpp::export]]
// CharacterVector cpp_par_newseqs(const IntegerMatrix freq_len, const CharacterVector seqs, 
//                                     const IntegerMatrix freq_mat, int n_samps) {
//     
//     // allocate the matrix we will return
//     CharacterVector output(seqs.size() * n_samps);
//     
//     // create the worker
//     new_seq new_seq(freq_len, seqs, freq_mat, n_samps, output);
//     
//     // call it with parallelFor
//     parallelFor(0, freq_len.nrow(), new_seq);
//     
//     return output;
// }








