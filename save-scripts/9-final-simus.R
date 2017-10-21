library(tidyverse)

n <- nrow(G)
n.train <- 6e3

params.grid <- expand.grid(
  par.causal = list(c(30, "all"), c(300, "all"), c(3000, "all"), c(30, "HLA")), 
  par.dist   = c("gaussian", "laplace"), 
  par.h2     = c(0.5, 0.8), 
  par.model  = c("simple", "fancy"),
  num.simu   = 1:100,
  stringsAsFactors = FALSE
) 



for (i in rows_along(results)) {
  
  res.file <- paste("results1/simu_", i)
  if (file.exists(res.file)) next
  
  params <- params.grid[i, ]
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
  ind.train <- sort(sample(n, size = N.train))
  ind.test <- setdiff(1:n, ind.train)
    
  # Get results from all methods
  res <- bind_rows(
    PRS(G, CHR, POS, pheno.all, covar.all, ind.train, ind.test),
    ttrees(TTree, "ttrees", pheno.all, ind.train, ind.test, n.trees = 100),
    logit.CMSA(G,  pheno.all, covar.all, ind.train, ind.test, "logit-simple"),
    logit.CMSA(G2, pheno.all, covar.all, ind.train, ind.test, "logit-triple")
  )
  print(res)
  saveRDS(res, file = res.file)
}