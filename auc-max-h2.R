nsimu <- 1e7
K <- 0.3

replicate(10, {
  sapply(c(0.5, 0.8), function(h2) {
    S <- rnorm(nsimu, mean = 0, sd = sqrt(h2))
    E <- rnorm(nsimu, mean = 0, sd = sqrt(1 - h2))
    Y <- (S + E) > qnorm(1 - K)
    
    bigstatsr::AUC(S, Y)
  })
})

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

