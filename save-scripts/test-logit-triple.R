system.time(
  G2 <- snp_triple(celiac2)
)


print(system.time(
  cmsa.logit3 <- big_CMSA(FUN = big_spLogReg, feval = AUC,
                          X = G2, y.train = pheno[ind.train], 
                          ind.train = ind.train, 
                          covar.train = obj.svd2$u[ind.train, ],
                          alpha = 0.5, dfmax = 20e3, ncores = NCORES)
)) 
# 91.7%, 19 min with 6 cores


preds2 <- predict(cmsa.logit3, X = G2, ind.row = ind.test, 
                  covar.row = obj.svd2$u[ind.test, ])
AUCBoot(preds2, pheno[ind.test])
AUC(preds2, pheno[ind.test])


test <- big_spLogReg(G2, y01.train = pheno[ind.train], 
             ind.train = ind.train, 
             covar.train = obj.svd2$u[ind.train, ],
             alpha = 0.5, dfmax = 20e3)

library(Matrix)
preds3 <- predict(test, X = G2, ind.row = ind.test, 
                  covar.row = obj.svd2$u[ind.test, ])
AUCBoot(preds3, pheno[ind.test])
AUC(preds3, pheno[ind.test])


counts <- big_counts(G2)
tmp <- big_scale()(G2)
