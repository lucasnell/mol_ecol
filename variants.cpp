#include <Rcpp.h>
#include <algorithm>

/*
 
 // RcppArmadillo version of sample function and its dependencies (deprecated for now)
 
 #include <RcppArmadilloExtensions/sample.h>
 // [[Rcpp::depends(RcppArmadillo)]]
 
 IntegerVector csample_int(IntegerVector x, int size,
                           bool replace = false, 
                           NumericVector prob = NumericVector::create()) {
     IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
     return ret;
 }
*/

using namespace Rcpp;
using namespace std;




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
 Random number generator
 (From here: http://xoroshiro.di.unimi.it/splitmix64.c)
 ------------
*/

uint32_t rng(void) {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;
    t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}

/*
 ------------
 Random integer generation, within a range
 ------------
*/

float rng_max = 2147483648; // 2^31

// [[Rcpp::export]]
IntegerVector int_sampler(int num_samps, float range_min, float range_max) {
    IntegerVector out_vec(num_samps);
    IntegerVector tmp_vec(num_samps);
    NumericVector num_vec;
    int x;
    for (int i=0; i < num_samps; i++) {
        x = rng() >> 1;
        tmp_vec[i] = x;
    }
    num_vec = tmp_vec;
    num_vec = (num_vec * (range_max - range_min) / rng_max) + range_min;
    out_vec = ceil(num_vec);
    return out_vec;
}




/*
 ------------
 Get the site locations for multiple sequences
 ------------
*/

// [[Rcpp::export]]
IntegerVector cpp_get_sites(IntegerVector seq_lens, IntegerMatrix seq_freq) {
    
    int n = seq_freq.nrow();
    int nn = sum(seq_freq(_,1));
    int seq_num;
    int samp_n;
    int seq_length;
    int j=0;
    
    IntegerVector ran_locs(nn);
    IntegerVector tmp_locs;
    
    for(int i = 0; i < n; i++) {
        seq_num = seq_freq(i,0) - 1;
        samp_n = seq_freq(i,1);
        seq_length = seq_lens[seq_num];
        tmp_locs = int_sampler(samp_n, 1, seq_length);
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
CharacterVector cpp_change_sites_1s(string seq, IntegerVector positions,
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
    ran_c = shuffle_int(ran_c);
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




// // Change sites for multiple sequences
// // [[Rcpp::export]]
// CharacterVector cpp_change_sites(CharacterVector seq_vec, IntegerMatrix positions_mat,
//                                  IntegerMatrix freq_mat) {
//     int n_samps = sum(freq_mat(0,_));
//     int n_seqs = seq_vec.size();
//     int out_seqs = n_seqs * n_samps;
//     CharacterVector out_vec(out_seqs);
//     string seq_i;
//     IntegerVector positions_i;
//     IntegerVector seq_nums = positions_mat(_,0);
//     IntegerVector positions = positions_mat(_,1);
//     CharacterVector new_seqs;
//     int j = 0;
// 
//     for (int i = 0; i < n_seqs; i++) {
//         seq_i = seq_vec[i];
//         positions_i = positions[seq_nums == (i + 1)];
//         new_seqs = cpp_change_sites_1s(seq_i, positions_i - 1, freq_mat, n_samps);
//         out_vec[Range(j, j + n_samps - 1)] = new_seqs;
//         j += n_samps;
//     }
// 
//     return out_vec;
// }
// 
