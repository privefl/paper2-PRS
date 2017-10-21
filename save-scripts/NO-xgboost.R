G2 <- celiac2$genotypes[, celiac2$map$chromosome == 22]

for (nr in c(10, 20, 50, 100)) {
  tmp <- xgboost::xgboost(data = G2[ind.train, ], label = pheno[ind.train],
                          nrounds = nr, verbose = 0)
  preds <- predict(tmp, newdata = G2[ind.test, ])
  print(AUCBoot(preds, pheno[ind.test]))
}


