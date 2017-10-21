tmp1 <- PRS(G, CHR, POS, pheno.all, covar.all, ind.train, ind.test)
map_dbl(tmp1$eval, ~AUC(.x[, 1], .x[, 2]))

tmp2 <- logit.CMSA(G,  pheno.all, covar.all, ind.train, ind.test, "logit-simple")


tmp3 <- logit.CMSA(G2, pheno.all, covar.all, ind.train, ind.test, "logit-triple")


tmp4 <- ttrees(TTree, "backingfiles/ttrees", pheno.all, ind.train, ind.test, n.trees = 10)

tmp <- bind_rows(tmp1, tmp2, tmp3, tmp4)
tmp$AUC <- map_dbl(tmp$eval, ~AUC(.x[, 1], .x[, 2]))
tmp

param1 <- params.grid[1, ]
param1$res <- list(tmp)
test <- unnest(param1, res)
bind_cols(param1[rep(1, 6), ], tmp)
