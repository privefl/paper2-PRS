/******************************************************************************/

// [[Rcpp::depends(RcppArmadillo, bigstatsr, BH)]]
#include <RcppArmadillo.h>
#include <triple-accessor.h>
#include <bigstatsr/biglasso/logistic.hpp>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
List triple_cdfit_binomial_hsr(Environment BM,
                               const NumericVector& y,
                               const IntegerVector& row_idx,
                               const NumericMatrix& covar,
                               NumericVector& lambda,
                               int L,
                               bool lam_scale,
                               double lambda_min,
                               double alpha,
                               bool user,
                               double eps,
                               int max_iter,
                               const NumericVector& m,
                               int dfmax,
                               bool warn,
                               bool verbose) {
  
  XPtr<FBM> xpBM = BM["address"];
  IntegerVector rows = row_idx - 1;
  
  return bigstatsr::biglassoLog::COPY_cdfit_binomial_hsr(
    RawSubTripleAcc(xpBM, rows, covar, BM["code256"]), 
    y, lambda, L, lam_scale, lambda_min, alpha,
    user, eps, max_iter, m, dfmax, warn, verbose
  );    
}

/******************************************************************************/