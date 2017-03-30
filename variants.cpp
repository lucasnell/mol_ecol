#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

IntegerVector csample_int(IntegerVector x, int size,
                          bool replace, 
                          NumericVector prob = NumericVector::create()) {
    IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
    return ret;
}


IntegerVector cpp_one_seq(int samp_n, int seq_length) {
    
    NumericVector ran_floats;
    int k=1;
    
    IntegerVector ran_locs(samp_n);
    
    ran_floats = runif(samp_n, 0, seq_length);
    ran_locs = ceiling(ran_floats);
    
    while (is_true(any(duplicated(ran_locs))) && k < 1000) {
        ran_floats = runif(samp_n, 0, seq_length);
        ran_locs = ceiling(ran_floats);
        k += 1;
    }
    return ran_locs;
}

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
        tmp_locs = cpp_one_seq(samp_n, seq_length);
        // tmp_locs = csample_int(Range(1, seq_length), samp_n, FALSE);
        ran_locs[Range(j, j + samp_n - 1)] = tmp_locs;
        j += samp_n;
    }
    
    return ran_locs;
}


