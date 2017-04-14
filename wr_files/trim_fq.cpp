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
int cpp_trim_fq(string in_fn, string out_fn) {
    string line, tmp_str;
    int line_num = 1;
    ifstream in_file(in_fn);
    ofstream out_file(out_fn);
    size_t ns;
    
    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            if (line_num == 2) {
                ns = line.find("N");
                if (ns != string::npos) {
                    tmp_str = line.substr(0, ns) + "\n";
                } else {
                    ns = line.size();
                    tmp_str = line + "\n";
                }
            } else if (line_num == 4) {
                tmp_str = line.substr(0, ns) + "\n";
                line_num = 0;
            } else {
                tmp_str = line + "\n";
            }
            out_file << tmp_str;
            line_num += 1;
        }
        in_file.close();
        out_file.close();
    } else {
        cout << "Unable to open file";
    }
    
    return 0;
}



// [[Rcpp::export]]
int cpp_trim_fqgz(string in_fn, string out_fn) {
    string line, tmp_str;
    int line_num = 1;
    size_t ns;
    // char* dptr;
    
    ifstream in_file(in_fn.c_str(), ios_base::in | ios_base::binary);
    ofstream out_file(out_fn);
    
    try {
        filtering_istream in_stream;
        in_stream.push(gzip_decompressor());
        in_stream.push(in_file);
        
        for(string line; getline(in_stream, line); ) {
            if (line_num == 2) {
                ns = line.find("N");
                if (ns != string::npos) {
                    tmp_str = line.substr(0, ns) + "\n";
                } else {
                    ns = line.size();
                    tmp_str = line + "\n";
                }
            } else if (line_num == 4) {
                tmp_str = line.substr(0, ns) + "\n";
                line_num = 0;
            } else {
                tmp_str = line + "\n";
            }
            out_file << tmp_str;
            line_num += 1;
        }
        in_file.close();
        out_file.close();
    }
    catch(const gzip_error& e) {
        cout << e.what() << '\n';
    }
    
    return 0;

}