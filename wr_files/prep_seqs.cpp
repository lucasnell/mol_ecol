#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

string rcomp(string seq) {
    
    string out_seq = seq;
    
    reverse(out_seq.begin(), out_seq.end());
    
    for (size_t i = 0; i < out_seq.length(); ++i){
        switch (out_seq[i]) {
        case 'A':
            out_seq[i] = 'T';
            break;
        case 'C':
            out_seq[i] = 'G';
            break;
        case 'G':
            out_seq[i] = 'C';
            break;
        case 'T':
            out_seq[i] = 'A';
            break;
        }
    }
    
    return out_seq;
}



// [[Rcpp::export]]
CharacterVector all_rm(CharacterVector all_s, int cut_len, string f_bc, string r_bc) {
    std::string f_strand, r_strand, s;
    int seq_len, vec_len = all_s.length();
    CharacterVector out_vec(vec_len * 2);
    for (int i = 0; i < vec_len; i++) {
        s = all_s[i];
        seq_len = s.size();
        if (seq_len < cut_len) {
            f_strand = s;
            r_strand = rcomp(s);
        } else {
            f_strand = s.substr(0, cut_len);
            r_strand = s.substr(seq_len - cut_len, cut_len);
            r_strand = rcomp(r_strand);
        }
        out_vec[i] = f_bc + f_strand;
        out_vec[vec_len + i] = r_bc + r_strand;
    }
    return out_vec;
}