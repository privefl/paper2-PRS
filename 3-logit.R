print(system.time(
  cmsa.logit <- big_CMSA(FUN = big_spLogReg, feval = AUC,
                         X = G, y.train = pheno[ind.train], 
                         ind.train = ind.train, 
                         covar.train = obj.svd2$u[ind.train, ],
                         alpha = 0.5, dfmax = 20e3,
                         K = 12, ncores = NCORES)
)) # 10.5 min with 6 core
