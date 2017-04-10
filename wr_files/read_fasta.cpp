#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace Rcpp;
using namespace std;
using namespace boost::iostreams;
// [[Rcpp::depends(BH)]]





// [[Rcpp::export]]
CharacterMatrix cpp_read_fasta(string file_name) {
    string line, tmp_name, tmp_seq;
    int line_num = 1;
    ifstream myfile(file_name);
    vector<string> seq_vec, name_vec;
    if (myfile.is_open()) {
        while (getline(myfile,line)) {
            if (line.find(">") != string::npos) {
                if (line_num > 1) {
                    seq_vec.push_back(tmp_seq);
                    name_vec.push_back(tmp_name);
                }
                line_num = 2;
                tmp_name = line.substr(1, line.size());
                tmp_seq = "";
            } else if (line != "\n") {
                tmp_seq += line.substr(0, line.size());
            }
        }
        myfile.close();
    } else {
        cout << "Unable to open file";
    }
    
    if (tmp_seq != "") {
        seq_vec.push_back(tmp_seq);
        name_vec.push_back(tmp_name);
    }
    
    CharacterMatrix out_mat(seq_vec.size(), 2);
    
    for (int i = 0; i < seq_vec.size(); i++) {
        out_mat(i,0) = name_vec[i];
        out_mat(i,1) = seq_vec[i];
    }
    
    return out_mat;
}


// [[Rcpp::export]]
CharacterMatrix cpp_read_fastagz(string file_name) {
    string line, tmp_name, tmp_seq;
    int line_num = 1;
    vector<string> seq_vec, name_vec;
    
    ifstream file(file_name.c_str(), ios_base::in | ios_base::binary);
    try {
        filtering_istream in_stream;
        in_stream.push(gzip_decompressor());
        in_stream.push(file);
        for(string line; getline(in_stream, line); ) {
            if (line.find(">") != string::npos) {
                if (line_num > 1) {
                    seq_vec.push_back(tmp_seq);
                    name_vec.push_back(tmp_name);
                }
                line_num = 2;
                tmp_name = line.substr(1, line.size());
                tmp_seq = "";
            } else if (line != "\n") {
                tmp_seq += line.substr(0, line.size());
            }
        }
        file.close();
    }
    catch(const gzip_error& e) {
        cout << e.what() << '\n';
    }

    if (tmp_seq != "") {
        seq_vec.push_back(tmp_seq);
        name_vec.push_back(tmp_name);
    }

    CharacterMatrix out_mat(seq_vec.size(), 2);

    for (int i = 0; i < seq_vec.size(); i++) {
        out_mat(i,0) = name_vec[i];
        out_mat(i,1) = seq_vec[i];
    }
    
    return out_mat;
    
}
