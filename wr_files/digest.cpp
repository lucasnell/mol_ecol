#include <Rcpp.h>
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
