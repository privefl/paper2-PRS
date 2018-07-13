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

/***R
add_sens_FDP3 <- function(x, corr, which_cond) {
  
  bigstatsr:::assert_class(corr, "dsCMatrix")
  corr <- as(abs(corr), "dgTMatrix")
  
  alt <- map2(x$true_set, x$set, ~ {
    if (length(.y) == 0) return(c(0, 0))
    bool.x <- bool.y <- logical(nrow(corr))
    bool.x[.x] <- TRUE
    bool.y[.y] <- TRUE
    ind <- which_cond(corr@i, corr@j, bool.x, bool.y)
    AltSens <-     sum(tapply(corr@x[ind], corr@i[ind], max)) / length(.x)
    AltFDP  <- 1 - sum(tapply(corr@x[ind], corr@j[ind], max)) / length(.y)
    c(AltSens, AltFDP)
  }) %>%
    transpose()
  
  x %>%
    mutate(AltSens = unlist(alt[[1]]),
           AltFDP  = unlist(alt[[2]]))
}
*/
