/******************************************************************************/

// [[Rcpp::depends(bigstatsr, BH)]]
#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
void tripleBM(Environment BM, Environment BM2) {
  
  XPtr<FBM> xpMat = BM["address"];
  int n = xpMat->nrow();
  int m = xpMat->ncol();
  SubBMCode256Acc macc(xpMat, seq_len(n)-1, seq_len(m)-1, BM["code256"]);
  
  XPtr<FBM> xpMat2 = BM2["address"];
  BMAcc<unsigned char> macc2(xpMat2);
  
  int i, j, j2;
  
  for (j = j2 = 0; j < m; j++, j2 += 3) {
    for (i = 0; i < n; i++) {
      macc2(i, j2)   = macc(i, j);
      macc2(i, j2+1) = macc(i, j) >= 0.5;
      macc2(i, j2+2) = macc(i, j) >  1.5;
    }
  }
}


/*** R
snp_triple <- function(x) {
  
  G0 <- x$genotypes
  G <- G0$copy(code = round(G0$code256))
  
  G2 <- FBM.code256(nrow(X), 3 * ncol(X), 
                    code = bigsnpr:::CODE_012,
                    backingfile = bigsnpr:::getNewFile(x, "tripled"))
  print(class(G2))
  tripleBM(G, G2)
  
  G2
}
*/
