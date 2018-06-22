library(tidyverse)
library(glue)
library(bigsnpr)
# https://stackoverflow.com/q/9486952/6103040
options(bigstatsr.cluster.type = "PSOCK")

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
  y.simu <- y.simu + rnorm(nrow(G), sd = sqrt(1 - h2))
  pheno <- as.numeric(y.simu > qnorm(1 - K))
}

celiac2 <- snp_attach("backingfiles/celiacQC_sub1.rds")
CHR <- celiac2$map$chromosome
POS <- celiac2$map$physical.pos
(G <- celiac2$genotypes)

covar.all <- readRDS("backingfiles/PCA2.rds")$u

n <- nrow(G)
ind.HLA <- snp_indLRLDR(CHR, POS, subset(LD.wiki34, ID == "hild12"))

NCORES <- nb_cores()


par.causal <- c(3000, "all")  ## 300 -> 3000
# Simulate phenotypes
pheno.all <- get_pheno(
  G,    
  h2 = 0.8, 
  M = as.integer(par.causal[1]), 
  ind.possible = `if`(par.causal[2] == "all", cols_along(G), ind.HLA),
  effects.dist = "gaussian", 
  model = "simple", 
  K = 0.3                 
)

# Split in training/test sets
ind.train <- sort(sample(n, size = 6000))
ind.test <- setdiff(1:n, ind.train)


system.time(
  cmsa.logit <- big_spLogReg(X = G, y01.train = pheno.all[ind.train], 
                             ind.train = ind.train, 
                             covar.train = covar.all[ind.train, , drop = FALSE],
                             alpha = 1, ncores = NCORES)
) # 44 sec -> 66 sec

preds.all <- rowMeans(
  predict(cmsa.logit, X = G, covar.row = covar.all, proba = FALSE)
)
AUC(preds.all[ind.test], pheno.all[ind.test])  ## 82.9% -> 57.5%


system.time(
  cmsa.logit2 <- big_spLogReg(X = G, y01.train = pheno.all[ind.train], 
                              base.train = preds.all[ind.train],
                              ind.train = ind.train, 
                              covar.train = covar.all[ind.train, , drop = FALSE],
                              alpha = 0.1, ncores = NCORES)
) # 20 sec -> 70 sec
preds.all2 <- rowMeans(
  predict(cmsa.logit2, X = G, covar.row = covar.all, proba = FALSE)
)
AUC(preds.all2[ind.test], pheno.all[ind.test]) # 0.5

system.time(
  cmsa.logit3 <- big_spLogReg(X = G, y01.train = pheno.all[ind.train], 
                              base.train = preds.all[ind.train],
                              ind.train = ind.train, 
                              covar.train = covar.all[ind.train, , drop = FALSE],
                              alpha = 0.001, ncores = NCORES)
) # 19 sec -> 97 sec
preds.all3 <- rowMeans(
  predict(cmsa.logit3, X = G, covar.row = covar.all, proba = FALSE)
)
AUC(preds.all3[ind.test], pheno.all[ind.test]) # 0.5

system.time(
  cmsa.logit4 <- big_spLogReg(X = G, y01.train = pheno.all[ind.train], 
                              # base.train = preds.all[ind.train],
                              ind.train = ind.train,
                              covar.train = covar.all[ind.train, , drop = FALSE],
                              alpha = 0.01, ncores = NCORES)
) # 102 sec
preds.all4 <- rowMeans(
  predict(cmsa.logit4, X = G, covar.row = covar.all, proba = FALSE)
)
AUC(preds.all4[ind.test], pheno.all[ind.test]) # 59.1%

system.time(
  cmsa.logit5 <- big_spLogReg(X = G, y01.train = pheno.all[ind.train], 
                              # base.train = preds.all[ind.train],
                              ind.train = ind.train,
                              covar.train = covar.all[ind.train, , drop = FALSE],
                              alpha = c(1, 0.001, 0.01, 1e-4), ncores = NCORES)
) # 218 sec
str(cmsa.logit5)
preds.all5 <- rowMeans(
  predict(cmsa.logit5, X = G, covar.row = covar.all, proba = FALSE)
)
AUC(preds.all5[ind.test], pheno.all[ind.test]) # 59.1%
