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

# Train logit-simple
system.time(
  logit <- big_spLogReg(X = G, y01.train = y[ind.train], 
                        ind.train = ind.train, 
                        covar.train = obj.svd$u[ind.train, ],
                        ncores = NCORES, alpha = c(1))
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


# Get predictions/individuals
preds1 <- predict(logit, X = G, ind.row = ind.test, 
                  covar.row = obj.svd$u[ind.test, ])
AUC(preds1, y[ind.test])

(nb.pred1 <- sum(rowSums(sapply(logit, function(l) l$beta.X != 0)) > 0))


ind.sets <- sample(rep_len(1:10, length(ind.train)))
test <- sapply(1:10, function(set) {
  ind.train <- ind.train[ind.sets != set]
  big_spLogReg(X = G, y01.train = y[ind.train], 
               ind.train = ind.train, 
               covar.train = obj.svd$u[ind.train, ],
               ncores = NCORES)
})
str(test)

ntimes <- rowSums(sapply(test, function(x) x$beta.X) > 0)
hist(ntimes)
table(ntimes)

test2 <- lapply(0:20, function(thr) {
  big_spLogReg(X = G, y01.train = y[ind.train], 
               ind.train = ind.train, 
               ind.col = which(ntimes >= thr),
               covar.train = obj.svd$u[ind.train, ],
               ind.sets = ind.sets,
               ncores = NCORES)
})
str(test2)

preds <- sapply(test2, function(mod) {
  predict(mod, G, ind.row = ind.test, covar.row = obj.svd$u[ind.test, ])
})
plot(apply(preds, 2, AUC, target = y[ind.test]), pch = 20)


sapply(test2, function(mod) {
  mean(rowSums(sapply(mod, function(x) x$beta.X) > 0) > 0)
})  ## still overfit
