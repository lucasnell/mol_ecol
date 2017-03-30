#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
IntegerVector cpp_one_seq(int i, NumericVector seq_lens, IntegerMatrix freq_mat) {
    int seq_length;
    NumericVector rand_seq_filt;
    int samp_n;
    // IntegerVector seq_range;
    int seq_num;
    NumericVector ran_floats;
    
    seq_num = freq_mat(i,0);
    samp_n = freq_mat(i,1);
    
    IntegerVector ran_locs(samp_n);

    seq_length = seq_lens[seq_num - 1];
    // seq_range = Range(1, seq_length);
    // ran_locs = RcppArmadillo::sample(seq_range, samp_n, FALSE);
    while (is_true(any(duplicated(ran_locs)))) {
        ran_floats = runif(samp_n, 0, seq_length);
        ran_locs = ceiling(ran_floats);
    }
    return ran_locs;
}

// IntegerVector cpp_get_sites(NumericVector seq_lens, NumericVector rand_seqs) {
//     int n = seq_lens.size();
//     
//     NumericVector rand_seq_filt;
//     IntegerVector ran_locs;
//     IntegerVector seq_range;
//     
//     for(int i = 0; i < n; i++) {
//         
//     }
//     
//     // seq_length = seq_lens[seq_num];
//     // seq_range = Range(1, seq_length);
//     // rand_seq_filt = rand_seqs[rand_seqs == (seq_num + 1)];
//     // samp_n = rand_seq_filt.size();
//     // ran_locs = RcppArmadillo::sample(seq_range, samp_n, FALSE);
//     // return ran_locs;
// }