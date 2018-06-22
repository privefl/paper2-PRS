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
                        ncores = NCORES, alpha = 0.001)
) 
## Alpha = 0.5
# FORK:  69 / 40 / 74
# PSOCK: 60 / 54 / 63
## Alpha = 1
# PSOCK: 58 / 65 / 97
## Alpha = 0.1
# PSOCK: 74 / 55 / 73
## Alpha = 0.01
# PSOCK: 1487 /  / 


# Get K = 10 predictions/individuals
preds1 <- predict(logit, X = G, ind.row = ind.test, 
                  covar.row = obj.svd$u[ind.test, ], proba = FALSE)
(nb.pred1 <- sum(rowSums(sapply(logit, function(l) l$beta.X != 0)) > 0))
(aucs1 <- apply(preds1, 2, AUC, target = y[ind.test]))
# Combining all K predictions
preds1.comb <- rowMeans(preds1)
AUC(preds1.comb, y[ind.test])
