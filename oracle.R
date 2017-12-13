get_oracle <- function(
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
  
  set <- sample(ind.possible, size = M)
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
  stopifnot(all.equal(drop(var(y.simu)), h2))
  S <- y.simu
  y.simu <- S + rnorm(nrow(G), sd = sqrt(1 - h2))
  pheno <- as.numeric(y.simu > qnorm(1 - K))
  AUC(S, pheno)
}

get_oracle(G, 0.8, 30, effects.dist = "gaussian", model = "simple")
get_oracle(G, 0.8, 3000, effects.dist = "gaussian", model = "simple")


params.grid2 <- expand.grid(
  par.causal = list(c(30, "all"), c(300, "all"), c(3000, "all"), c(30, "HLA")), 
  par.dist   = c("gaussian", "laplace"), 
  par.h2     = c(0.5, 0.8), 
  par.model  = c("simple", "fancy"),
  num.simu   = 1:100,
  stringsAsFactors = FALSE
) 

AUC_oracle <- function(params) {
  
  cat(".")
  
  par.causal <- params[["par.causal"]][[1]]
  
  # Simulate phenotypes
  AUC <- get_oracle(
    G,    
    h2 = params[["par.h2"]], 
    M = as.integer(par.causal[1]), 
    ind.possible = `if`(par.causal[2] == "all", cols_along(G), ind.HLA),
    effects.dist = params[["par.dist"]], 
    model = params[["par.model"]], 
    K = 0.3                 
  )
  
  params$AUC <- AUC
  params
}
oracles <- lapply(rows_along(params.grid2), function(i) AUC_oracle(params.grid2[i, ]))
oracles2 <- do.call(bind_rows, oracles)
oracles2 %>% 
  mutate(
    par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all")))) %>%
  group_by(par.causal, par.dist, par.h2, par.model) %>%
  summarise(AUC_mean = round(mean(AUC), 3),
            AUC_mean_boot = boot(AUC, 1e4, mean))
