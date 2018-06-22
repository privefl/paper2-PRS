# https://www.ncbi.nlm.nih.gov/pubmed/27219331
for (N in c(50, 500, 5000)) {
  seq_P <- exp(seq(log(50), log(5000), length.out = 20))
  res <- sapply(seq_P, function(p) {
    X <- matrix(rnorm(N * P), N)
    corr <- abs(cor(X[, 1:5], X[, -(1:5)]))
    
    AltSens <- mean(apply(corr, 1, max))
    AltFDP <- 1 - mean(apply(corr, 2, max))
    c(AltSens, AltFDP)
  })
  
  plot(seq_P, res[1, ], pch = 20, ylim = c(0, 1), log = 'x')
  points(seq_P, res[2, ], pch = 20, col = 3)
}


