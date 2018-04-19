params.grid2 <- expand.grid(
  n.train = 3000,
  par.causal = list(c(30, "all"), c(300, "all"), c(3000, "all"), c(30, "HLA")), 
  par.dist   = c("gaussian", "laplace"), 
  par.h2     = c(0.5, 0.8), 
  par.model  = c("simple", "fancy"),
  num.simu   = 1:100,
  stringsAsFactors = FALSE
) 

i <- 1
params <- params.grid2[i, ]
par.causal <- params[["par.causal"]][[1]]

# Simulate phenotypes
pheno.all <- get_pheno(
  G,    
  h2 = params[["par.h2"]], 
  M = as.integer(par.causal[1]), 
  ind.possible = `if`(par.causal[2] == "all", cols_along(G), ind.HLA),
  effects.dist = params[["par.dist"]], 
  model = params[["par.model"]], 
  K = 0.3                 
)

# Split in training/test sets
ind.train <- sort(sample(n, size = params[["n.train"]]))
ind.test <- setdiff(1:n, ind.train)

system.time(
  test <- logit.CMSA(G, pheno.all, covar.all, 
                     ind.train, ind.test, "logit-simple")
)
AUC(test$eval[[1]][, 1], test$eval[[1]][, 2])

library(biglasso)
# G2 <- big.matrix(nrow(G), ncol(G) + ncol(covar.all), 
#                 backingfile = "G-PC", backingpath = "backingfiles")
# dim(G2)
# big_apply(G, function(X, ind) {
#   G2[, ind] <- X[, ind]
#   NULL
# }, a.combine = 'c')
# G2[, ncol(G) + cols_along(covar.all)] <- covar.all
G2 <- attach.big.matrix("backingfiles/G-PC.desc")
system.time(
  test2 <- biglasso(G2, pheno.all, ind.train, ncores = nb_cores(),
                    penalty = "enet", family = "binomial", alpha = 0.5)
)

