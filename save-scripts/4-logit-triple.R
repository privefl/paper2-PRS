
print(system.time(
  cmsa.logit3 <- big_CMSA(FUN = triple_biglasso, feval = AUC,
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

