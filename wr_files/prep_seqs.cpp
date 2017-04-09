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
CharacterVector prep_s(CharacterVector all_s, int f_len, int r_len, 
                       string f_bc, string r_bc) {
    std::string f_strand, r_strand, s;
    int frag_len, vec_len = all_s.length();
    CharacterVector out_vec(vec_len * 2);
    for (int i = 0; i < vec_len; i++) {
        s = all_s[i];
        frag_len = s.size();
        if (frag_len >= r_len & frag_len >= f_len) {
            f_strand = s.substr(0, f_len);
            r_strand = s.substr(frag_len - r_len, r_len);
            r_strand = rcomp(r_strand);
        } else if (frag_len < r_len & frag_len < f_len) {
            f_strand = s;
            r_strand = rcomp(s);
        } else if (frag_len >= r_len) {
            f_strand = s;
            r_strand = s.substr(frag_len - r_len, r_len);
            r_strand = rcomp(r_strand);
        } else {
            f_strand = s.substr(0, f_len);
            r_strand = rcomp(s);
        }
        out_vec[i] = f_bc + f_strand;
        out_vec[vec_len + i] = r_bc + r_strand;
    }
    return out_vec;
}


// [[Rcpp::export]]
CharacterVector prep_p(CharacterVector all_s, int f_len, int r_len, 
                       string f_bc, string r_bc) {
    string f_strand, r_strand, s, insert;
    int frag_len;
    int insert_len = f_len + r_len; // insert length without barcodes
    int vec_len = all_s.length();
    CharacterVector out_vec(vec_len);
    for (int i = 0; i < vec_len; i++) {
        s = all_s[i];
        frag_len = s.size();
        if (frag_len > insert_len) {
            f_strand = s.substr(0, f_len);
            r_strand = s.substr(frag_len - r_len, r_len);
            insert = f_bc + f_strand + r_strand + r_bc;
        } else {
            insert = f_bc + s + r_bc;
        }
        out_vec[i] = insert;
    }
    return out_vec;
}