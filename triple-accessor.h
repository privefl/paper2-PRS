#ifndef SUBTRIPLEACC_H
#define SUBTRIPLEACC_H

/******************************************************************************/

#include <bigstatsr/SubMatCovAcc.h>

/******************************************************************************/

class RawSubTripleAcc : public RawSubMatCovAcc {
public:
  RawSubTripleAcc(const FBM * xpBM,
                  const IntegerVector& row_ind,
                  const NumericMatrix& covar,
                  const NumericVector& code256)
    : RawSubMatCovAcc(xpBM, row_ind, covar, code256) {
    
      NumericMatrix tmp(code256.size(), 3);
      tmp(_, 0) = code256;
      tmp(_, 1) = code256 >= 0.5;
      tmp(_, 2) = code256 > 1.5;
      _lookup = tmp;
    }
  
  inline double operator() (size_t i, size_t j) {
    if (j < ncol()) {
      int j2 = j / 3; 
      return _lookup(SubMatCovAcc<unsigned char>::operator()(i, j2), j % 3);
    } else {
      return _covar(i, j - ncol());
    }
  }
  
  size_t ncol() const {
    return 3 * _ncol + _ncoladd;
  }
  
protected:
  NumericMatrix _lookup;
};

/******************************************************************************/

#endif // SUBTRIPLEACC_H