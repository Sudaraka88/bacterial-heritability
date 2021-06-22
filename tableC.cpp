#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

struct char_cmp { 
  bool operator () (const char *a,const char *b) const 
  {
    return strcmp(a,b)<0;
  } 
};
typedef std::map<const char *, int, char_cmp> Map;
// [[Rcpp::export]]
std::map<const char *, int, char_cmp> tableC(StringVector x) {
  Map m;
  m.clear();
  // Fill the map with occurences per number
  for (int i = 0; i != x.size(); ++i) {
    m[x[i]] += 1;
  }
  return m;
}

// Hints from https://thecoatlessprofessor.com/programming/cpp/porting-rs-table-function-to-c-/
