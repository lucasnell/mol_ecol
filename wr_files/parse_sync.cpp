#include <Rcpp.h>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <algorithm>



using namespace Rcpp;
using namespace std;



/*
 Split a string by a delimiter
 
 "The first puts the results in a pre-constructed vector, the second returns a new 
 vector." (http://stackoverflow.com/a/236803/5016095)
 
 ** This only works for single-character delimiters **
 
 */
template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vector<string> cpp_split_delim(const string &s, char delim = '\t') {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

/*
 This returns a list containing nearly all the info present in the sync file
 _Note_: n_samps is # samples in addition to pooled sample
 */
// [[Rcpp::export]]
List cpp_parse_sync_all_info(string in_fn, int n_samps) {
    string line;
    vector<string> line_vec, tmp_dep, tmp_dep_row, scaffolds;
    int dep_cols = (n_samps + 1) * 4;
    int tmp_pos;
    vector<int> positions, dep_row(dep_cols);
    vector<vector<int> > depths;
    ifstream in_file(in_fn);

    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            line_vec = cpp_split_delim(line);
            scaffolds.push_back(line_vec[0]);
            stringstream(line_vec[1]) >> tmp_pos;
            positions.push_back(tmp_pos);
            tmp_dep_row = cpp_split_delim(line_vec[3], ':');
            tmp_dep_row.resize(4);
            for (int i = 4; i < (n_samps + 4); i++) {
                tmp_dep = cpp_split_delim(line_vec[i], ':');
                tmp_dep.resize(4);
                tmp_dep_row.insert(tmp_dep_row.end(), tmp_dep.begin(), tmp_dep.end());
            }
            for (int j = 0; j < dep_cols; j++) {
                stringstream(tmp_dep_row[j]) >> dep_row[j];
            }
            depths.push_back(dep_row);
        }
        in_file.close();
    } else {
        stop("Unable to open file");
    }
    
    IntegerMatrix dep_mat(depths.size(), dep_cols);
    IntegerVector tmp_vec;
    for(int i = 0; i < depths.size(); i ++) {
        tmp_vec = wrap(depths[i]);
        dep_mat(i, _) = tmp_vec;
    }
    
    List ret_list;
    ret_list["scaffolds"] = scaffolds;
    ret_list["positions"] = positions;
    ret_list["depths"] = dep_mat;
    return ret_list;
}



// This gets the positions of nucleotides that could have been contributed to the 
// pooled DNA sequence
// Input is a 4-length vector of frequencies of nucleotides A, T, C, G (in order)
vector<int> cpp_pool_pos(vector<int> pool_freqs, int min_value) {
    
    vector<int> out_pos;
    
    for (int i = 0; i < pool_freqs.size(); i++) {
        if (pool_freqs[i] >= min_value) {
            out_pos.push_back(i);
        }
    }
    return out_pos;
}

// This gets the number of alleles that could have contributed to the pooled DNA 
// sequence
// This assumes diploid
// Both inputs are 4-length vectors of frequencies of nucleotides A, T, C, G (in order)
int cpp_allele_contr(vector<int> pool_poss, vector<int> samp_freqs, int min_freq, 
                     double het_cutoff, bool diploid) {
    
    int alleles = 0;
    bool in_pool;
    bool hetero;
    
    // First get indices of first and second-highest frequencies
    int max_f = samp_freqs[0];
    int smax_f = 0;
    int max_i = 0;
    int smax_i = 0;
    for(int i = 1; i < samp_freqs.size(); i++) {
        if (samp_freqs[i] > max_f){
            smax_f = max_f;
            smax_i = max_i;
            max_f = samp_freqs[i];
            max_i = i;
        } else if (samp_freqs[i] > smax_f){
            smax_f = samp_freqs[i];
            smax_i = i;
        }
    }
    // If the maximum frequency isn't above cutoff, return 0
    if (max_f < min_freq) {
        return 0;
    }

    in_pool = find(pool_poss.begin(), pool_poss.end(), max_i) != pool_poss.end();
    // If the max index is found in the pool possible vector, assume the max-frequency
    // allele could be contributed to the pool
    if (in_pool) {
        alleles += 1;
    }
    
    // If the organism is diploid, the ways it can contribute two alleles are...
    // 1. It's heterozygous and the second-highest frequency allele is found in the pool
    // 2. It's homozygous and both alleles are found in the pool (which we've already
    //    checked for)
    if (diploid) {
        hetero = smax_f >= (max_f * het_cutoff);
        if (hetero) {
            in_pool = find(pool_poss.begin(), pool_poss.end(), smax_i) != pool_poss.end();
            if (in_pool) {
                alleles += 1;
            }
        } else if (in_pool) {
            alleles += 1;
        }
    }

    return alleles;
}



/*
 This returns a list containing much less info present than in the sync file
 _Note #1_: n_samps is # samples in addition to pooled sample
 _Note #2_: It is assumed that the first sample in the sync file is the pooled sample
*/
// [[Rcpp::export]]
List cpp_parse_sync_concise(string in_fn, int n_samps, int min_pool = 0, 
                            int min_freq = 10, double het_cutoff = 0.5, 
                            bool diploid = true) {
    // Reading lines from file
    string line;
    ifstream in_file(in_fn);
    vector<string> line_vec;
    // Pooled DNA read depth
    vector<string> tmp_pool_freqs;
    vector<int> pool_freqs(4);
    vector<vector<int> > pool_depths;
    // 2 vectors indicating position on genome: scaffolds and positions
    vector<string> scaffolds;
    int tmp_pos;
    vector<int> positions;
    // Getting possible sample allele contributions
    vector<int> samp_freqs(4);
    vector<int> pool_poss, samp_row(n_samps);
    vector<string> tmp_samp_freqs;
    vector<vector<int> > samp_alleles;

    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            // Split by tab
            line_vec = cpp_split_delim(line);
            // Append scaffold name
            scaffolds.push_back(line_vec[0]);
            // Append position as integer
            stringstream(line_vec[1]) >> tmp_pos;
            positions.push_back(tmp_pos);
            // Append pooled DNA nucleotide frequencies as integers
            tmp_pool_freqs = cpp_split_delim(line_vec[3], ':');
            for (int j = 0; j < 4; j++) {
                stringstream(tmp_pool_freqs[j]) >> pool_freqs[j];
            }
            pool_depths.push_back(pool_freqs);
            // Append sample nucleotides
            // First finding positions of nucleotides that were contributed to pool:
            pool_poss = cpp_pool_pos(pool_freqs, min_pool);
            if (pool_poss.size() == 0) {
                fill(samp_row.begin(), samp_row.end(), 0);
            // If >0 pooled nucleotide type, then fill in possible contributions for
            // each sample
            } else {
                for (int i = 0; i < n_samps; i++) {
                    // Getting frequencies like for the pooled sample
                    tmp_samp_freqs = cpp_split_delim(line_vec[(i + 4)], ':');
                    for (int j = 0; j < 4; j++) {
                        stringstream(tmp_samp_freqs[j]) >> samp_freqs[j];
                    }
                    // Calculating possible allele contributions
                    samp_row[i] = cpp_allele_contr(pool_poss, samp_freqs, min_freq, 
                                                   het_cutoff, diploid);
                }
            }
            samp_alleles.push_back(samp_row);
        }
        in_file.close();
    } else {
        stop("Unable to open file");
    }

    IntegerMatrix pool_mat(pool_depths.size(), 4);
    IntegerVector tmp_pvec;
    for(int i = 0; i < pool_depths.size(); i ++) {
        tmp_pvec = wrap(pool_depths[i]);
        pool_mat(i, _) = tmp_pvec;
    }

    NumericMatrix alleles_mat(samp_alleles.size(), n_samps);
    NumericVector tmp_avec;
    // int tmp_sum;
    vector<bool> all_equal_vec(samp_alleles.size());
    bool all_equal;
    for(int i = 0; i < samp_alleles.size(); i ++) {
        tmp_avec = wrap(samp_alleles[i]);
        // tmp_sum = sum(tmp_avec);
        // tmp_sum = std::max(1, tmp_sum);
        // tmp_avec = tmp_avec / tmp_sum;
        tmp_avec = tmp_avec;
        alleles_mat(i, _) = tmp_avec;
        all_equal = adjacent_find(tmp_avec.begin(), tmp_avec.end(), 
                                  not_equal_to<int>()) == tmp_avec.end();
        all_equal_vec[i] = all_equal;
    }
    List ret_list;
    ret_list["scaffolds"] = scaffolds;
    ret_list["positions"] = positions;
    ret_list["pool_depths"] = pool_mat;
    ret_list["samp_alleles"] = alleles_mat;
    ret_list["all_equal"] = all_equal_vec;
    return ret_list;
}



// [[Rcpp::export]]
List cpp_get_genome(string in_fn) {
    
    ifstream in_file(in_fn);
    string line;
    vector<string> names;
    vector<int> sizes;
    int tmp_size, line_num = 1;
    string tmp_name;
    
    
    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            if (line.find(">") != string::npos) {
                if (line_num > 1) {
                    names.push_back(tmp_name);
                    sizes.push_back(tmp_size);
                }
                line_num = 2;
                tmp_name = line.substr(1, line.find(' ') - 1);
                tmp_size = 0;
            } else if (line != "\n") {
                tmp_size += line.size();
            }
        }
        in_file.close();
    } else {
        stop("Unable to open file");
    }
    
    names.push_back(tmp_name);
    sizes.push_back(tmp_size);
    
    List ret_list;
    ret_list["names"] = names;
    ret_list["sizes"] = sizes;
    return ret_list;
}
