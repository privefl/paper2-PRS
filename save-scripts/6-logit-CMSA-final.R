logit.CMSA <- function(G, pheno.all, covar.all, ind.train, ind.test, method) {
  
  timing <- system.time({
    
    cmsa.logit <- big_CMSA(FUN = big_spLogReg, feval = AUC,
                           X = G, y.train = pheno.all[ind.train], 
                           ind.train = ind.train, 
                           covar.train = covar.all[ind.train, , drop = FALSE],
                           alpha = 0.5, dfmax = 20e3, 
                           ncores = nb_cores())
    
    preds <- predict(cmsa.logit, X = G, ind.row = ind.test, 
                     covar.row = covar.all[ind.test, , drop = FALSE])
  })[3]
  
  tibble(
    method   = method, 
    eval     = list(cbind(preds, pheno.all[ind.test])),
    timing   = timing,
    nb.preds = sum(cmsa.logit != 0)
  )
}

test <- logit.CMSA(G, pheno, obj.svd2$u, ind.train, ind.test, "LOGIT")
test2 <- logit.CMSA(G2, pheno, obj.svd2$u, ind.train, ind.test)

str(test)
AUC(test$eval[, 1], test$eval[, 2])
test$timing
test$infos

str(test2)
AUC(test2$eval[, 1], test2$eval[, 2])
test2$timing
test2$infos
