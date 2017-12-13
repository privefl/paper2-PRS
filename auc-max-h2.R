
K <- 0.3

nsimu <- 1e8
replicate(5, {
  sapply(c(0.5, 0.8), function(h2) {
    S <- rnorm(nsimu, mean = 0, sd = sqrt(h2))
    E <- rnorm(nsimu, mean = 0, sd = sqrt(1 - h2))
    Y <- (S + E) > qnorm(1 - K)
    
    bigstatsr::AUC(S, Y)
  })
})
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.8410003 0.8410194 0.8410863 0.8410095 0.8410279
# [2,] 0.9406166 0.9406530 0.9406010 0.9405916 0.9406338

# Estimation from Wray 2010
sapply(c(0.5, 0.8), function(h2) {
  T0 <- qnorm(1 - K)
  z <- dnorm(T0)
  i <- z / K
  v <- -i * K / (1 - K)
  num <- (i - v) * h2
  deno.part <- 1 - h2 * i * (i - T0) + (1 - h2 * v * (v - T0))
  pnorm(num / sqrt(h2 * deno.part))
})

