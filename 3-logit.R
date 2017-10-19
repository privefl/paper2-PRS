# Compute PCA
obj.svd2 <- snp_autoSVD(G, CHR, POS, ncores = NCORES)


print(system.time(
  cmsa.logit <- big_CMSA(FUN = big_spLogReg, feval = AUC,
                         X = G, y.train = pheno[ind.train], 
                         ind.train = ind.train, 
                         covar.train = obj.svd2$u[ind.train, ],
                         alpha = 0.5, dfmax = 20e3, ncores = NCORES)
)) 
# 91%, 7 min with 6 cores
# error, 31 min with 1 core


preds <- predict(cmsa.logit, X = G, ind.row = ind.test, 
                 covar.row = obj.svd2$u[ind.test, ])
AUCBoot(preds, pheno[ind.test])
AUC(preds, pheno[ind.test])


n <- length(ind.train)
K <- 10
indCV <- sample(rep_len(1:K, n))
cl <- parallel::makeCluster(NCORES)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)
library(foreach)
cross.res <- foreach(ic = 1:K) %dopar% {
  in.val <- (indCV == ic)
  mod <- FUN(X, y.train[!in.val], ind.train[!in.val], covar.train[!in.val,
                                                                  ],
             alpha = 0.5, dfmax = 20e3)
  scores <- predict(mod, X, ind.row = ind.train[in.val],
                    covar.row = covar.train[in.val, ])
  list(betas = mod$beta, scores = scores)
}
betas <- sapply(cross.res, function(x) {
  tmp <- x$scores
  ind <- as.numeric(rownames(tmp))
  stopifnot(all(ind %in% ind.train))
  ind2 <- match(ind, ind.train)
  seval <- apply(tmp, 2, feval, target = y.train[ind2])
  x$betas[, which.max(seval)]
})
beta <- bigstatsr:::get_beta(betas, "geometric-median")
ind <- which(beta != 0)
plot(ind, beta[ind])
