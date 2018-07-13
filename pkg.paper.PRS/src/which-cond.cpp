#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector which_cond(IntegerVector i, IntegerVector j,
                         LogicalVector x, LogicalVector y) {

  std::vector<int> res;
  int K = i.size();
  for (int k = 0; k < K; k++) {
    if (x[i[k]] && y[j[k]]) res.push_back(k + 1);
  }

  IntegerVector res2(res.size());
  std::copy(res.begin(), res.end(), res2.begin());
  return res2;
}
