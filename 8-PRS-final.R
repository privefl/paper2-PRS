PRS <- function(G, infos.chr, infos.pos, 
                pheno.all, covar.all, 
                ind.train, ind.test) {
  
  # GWAS
  gwas.train <- big_univLogReg(
    G, pheno.all[ind.train], ind.train = ind.train, 
    covar.train = covar.all[ind.train, , drop = FALSE], 
    ncores = nb_cores()
  )
  gwas.train.gc <- snp_gc(gwas.train)
  
  timing <- system.time({
    
    # Clumping on the test set
    ind.keep <- snp_clumping(G, infos.chr = infos.chr,
                             ind.row = ind.test,
                             thr.r2 = 0.2, 
                             S = abs(gwas.train.gc$score),
                             size = 500,
                             is.size.in.bp = TRUE,
                             infos.pos = infos.pos,
                             ncores = nb_cores())
    
    # PRS
    thrs <- c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100)))
    lpS <- -predict(gwas.train.gc)
    nb.pred <- sapply(thrs, function(thr) sum(lpS[ind.keep] > thr))
    prs <- snp_PRS(G, betas.keep = gwas.train.gc$estim[ind.keep],
                   ind.test = ind.test,
                   ind.keep = ind.keep,
                   lpS.keep = lpS[ind.keep], 
                   thr.list = thrs)
    
    preds <- predict(cmsa.logit, X = G, ind.row = ind.test, 
                     covar.row = covar.all[ind.test, , drop = FALSE])
  })[3]
  
  ind.best <- which.max(apply(prs, 2, AUC, target = pheno.all[ind.test]))
  
  methods <- c("PRS-all", "PRS-stringent", "PRS-max")
  indices <- c(1:2, ind.best)
  
  lapply(1:3, function(i) {
    k <- indices[i]
    tibble(
      method = methods[i],
      eval = list(cbind(prs[, k], pheno.all[ind.test])),
      timing = timing,
      nb.preds = nb.pred[k]
    )
  }) %>%
    bind_rows()
}
