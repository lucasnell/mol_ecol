#include <Rcpp.h>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>


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

// [[Rcpp::export]]
vector<string> cpp_split_delim(const string &s, char delim = '\t') {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}



// n_samps is # samples in addition to pooled sample
// [[Rcpp::export]]
List cpp_parse_sync(string in_fn, int n_samps) {
    string line, tmp_str;
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


