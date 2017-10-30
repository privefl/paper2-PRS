# Simulate some phenotypes
pheno.all <- get_pheno(G, h2 = 0.8, M = 30,
                       ind.possible = ind.HLA,
                       effects.dist = "gaussian", 
                       model = "simple")

# Split in training/test sets
ind.train <- sort(sample(n, size = n.train))
ind.test <- setdiff(1:n, ind.train)

# GWAS
gwas.train <- big_univLogReg(
  G, pheno.all[ind.train], ind.train = ind.train, 
  covar.train = covar.all[ind.train, , drop = FALSE], 
  ncores = nb_cores()
)
gwas.train.gc <- snp_gc(gwas.train)

# Clumping on the test set
ind.keep <- snp_clumping(G, infos.chr = CHR,
                         ind.row = ind.test,
                         thr.r2 = 0.2, 
                         S = abs(gwas.train.gc$score),
                         size = 500,
                         is.size.in.bp = TRUE,
                         infos.pos = POS,
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

aucs <- apply(prs, 2, AUC, target = pheno.all[ind.test])

cmsa.logit <- big_spLogReg(X = G, y01.train = pheno.all[ind.train], 
                           ind.train = ind.train, 
                           covar.train = covar.all[ind.train, , drop = FALSE],
                           alpha = 0.5, dfmax = 20e3)
preds <- predict(cmsa.logit, X = G, ind.row = ind.test, 
                 covar.row = covar.all[ind.test, , drop = FALSE])
aucs2 <- apply(preds, 2, AUC, target = pheno.all[ind.test])

library(Matrix)
bind_rows(
  data.frame(nb = nb.pred, AUC = aucs, method = "PRS"),
  data.frame(nb = colSums(cmsa.logit$beta != 0), AUC = aucs2, method = "logit")
) %>%
  myggplot(aes(nb, AUC, color = method)) +
  geom_point() + 
  geom_line() + 
  scale_x_log10()
