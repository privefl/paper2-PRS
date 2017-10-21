tmp1 <- PRS(G, CHR, POS, pheno.all, covar.all, ind.train, ind.test)
map_dbl(tmp1$eval, ~AUC(.x[, 1], .x[, 2]))

tmp2 <- logit.CMSA(G,  pheno.all, covar.all, ind.train, ind.test, "logit-simple")


tmp3 <- logit.CMSA(G2, pheno.all, covar.all, ind.train, ind.test, "logit-triple")


tmp4 <- ttrees(TTree, "ttrees", pheno.all, ind.train, ind.test, n.trees = 10)