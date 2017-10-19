get_pheno <- function(
  G,                                        ## matrix of genotypes
  h2,                                       ## heritability 
  M,                                        ## nbs of causal variants
  ind.possible = cols_along(G),             ## indices of possible causal variants
  effects.dist = c("gaussian", "laplace"),  ## distribution of effects 
  model = c("simple", "fancy"),             ## model for simulation
  K = 0.3                                   ## prevalence
) {
  
  effects.dist  <- match.arg(effects.dist)
  model <- match.arg(model)
  
  print(set <- sample(ind.possible, size = M))
  effects <- `if`(effects.dist == "gaussian", 
                  rnorm(M, sd = sqrt(h2 / M)),
                  rmutil::rlaplace(M, s = sqrt(h2 / (2*M))))
  
  if (model == "simple") {
    # only linear
    y.simu <- scale(G[, set]) %*% effects
  } else {
    sets <- split(1:M, sample(rep_len(1:3, M)))
    # linear
    ind1 <- sets[[1]]
    y.simu <- scale(G[, set[ind1]]) %*% effects[ind1]
    # recessive / dominant
    ind2 <- sets[[2]]
    y.simu <- y.simu + scale(G[, set[ind2]] > 0.5) %*% effects[ind2]
    # interactions
    ind3 <- matrix(sets[[3]], ncol = 2)
    y.simu <- y.simu + scale(G[, set[ind3[, 1]]] * G[, set[ind3[, 2]]]) %*% 
      effects[ind3[, 1]] * sqrt(2)
  }
  
  y.simu <- y.simu / sd(y.simu) * sqrt(h2)
  stopifnot(all.equal(drop(var(y.simu)), 0.8))
  y.simu <- y.simu + rnorm(nrow(G), sd = sqrt(1 - h2))
  pheno <- as.numeric(y.simu > qnorm(1 - K))
}

print(get_pheno(G, 0.8, 300, ind.possible = 
                  snp_indLRLDR(CHR, POS, subset(LD.wiki34, ID == "hild12"))))

##### Pas besoin de fonction pour ca
get_split <- function(n = nrow(G), n.train = 6e3) {
  ind.train <- sort(sample(n, size = N.train))
  ind.test <- setdiff(1:n, ind.train)
}
split(1:M, sample(rep_len(1:3, M)))

ind.possible <- snp_indLRLDR(CHR, POS, subset(LD.wiki34, ID == "hild12"))

