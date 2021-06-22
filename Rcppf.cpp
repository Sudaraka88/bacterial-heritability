#include <Rcpp.h>
#include <string>
#include <stdlib.h>
using namespace Rcpp;

// [[Rcpp::export]]
void makeFreqTable2(NumericVector N_ALLs, CharacterMatrix freq_only, int iter, NumericMatrix& freq_out) {
  // NumericMatrix freq_out(iter, 5);
  std::string temp = "";
  for(int i = 0; i < iter; i++){
    for(int j = 0; j < N_ALLs(i); j++){
      temp = as<std::string>(freq_only(i, j));
      std::string nuc = temp.substr(0,1);
      std::string freq = temp.substr(2,temp.length());
      // A C G T N order - flipped for faster performance!
      if(nuc == "C"){
        freq_out(i,1) = atof(freq.c_str());
      } else if(nuc == "G"){
        freq_out(i,2) = atof(freq.c_str());
      }else if(nuc == "T"){
        freq_out(i,3) = atof(freq.c_str());
      }else if(nuc == "A"){
        freq_out(i,0) = atof(freq.c_str());
      }else if(nuc == "*"){
        freq_out(i,4) = atof(freq.c_str());
        }
      }
    }
  // return(freq_out);
}
