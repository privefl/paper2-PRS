library(bigsnpr)
celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
NCORES <- nb_cores()
y <- celiac$fam$affection - 1
obj.svd <- readRDS("backingfiles/PCA.rds")

dim(G)
set.seed(1)
ind.train <- sort(sample(nrow(G), size = 12e3))
ind.test <- setdiff(rows_along(G), ind.train)

options(bigstatsr.cluster.type = "PSOCK")
# Train logit-simple
system.time(
  logit <- big_spLogReg(X = G, y01.train = y[ind.train], 
                        ind.train = ind.train, 
                        covar.train = obj.svd$u[ind.train, ],
                        ncores = NCORES, alpha = 1)
) # 41 sec
preds.all <- rowMeans(predict(logit, G, covar.row = obj.svd$u, proba = FALSE))
AUC(preds.all[ind.test], y[ind.test])

system.time(
  logit2 <- big_spLogReg(X = G, y01.train = y[ind.train], 
                         base.train = preds.all[ind.train],
                         ind.train = ind.train, 
                         covar.train = obj.svd$u[ind.train, ],
                         ncores = NCORES, alpha = 0.001)
) 
# 0.1   -> 21 sec
# 0.01  -> 21 sec
# 0.001 -> 21 sec


# Get K = 10 predictions/individuals
preds2 <- predict(logit2, X = G, ind.row = ind.test, 
                  covar.row = obj.svd$u[ind.test, ], proba = FALSE)
preds2 <- rowMeans(preds2)
AUC(preds2, y[ind.test])  # 0.5 à chaque fois
library(ggplot2)
qplot(preds2, fill = as.factor(y[ind.test]), geom = "density", alpha = I(0.4))

