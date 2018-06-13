library(bigstatsr)

# simulating some data
N <- 230
M <- 730
X <- FBM(N, M, init = rnorm(N * M, sd = 5))
y <- rowSums(X[, 1:5]) + 5 * rnorm(N)
covar <- matrix(rnorm(N * 3), N)

ind.train <- sort(sample(nrow(X), 150))
ind.test <- setdiff(rows_along(X), ind.train)

test <- big_spLinReg(X, y[ind.train], ind.train = ind.train,
                     covar.train = covar[ind.train, ],
                     alpha = 1)
# K = 10 predictions
str(preds <- predict(test, X, ind.row = ind.test, covar.row = covar[ind.test, ]))
# Combine them
preds2 <- rowMeans(preds)

plot(preds2, y[ind.test], pch = 20); abline(0, 1, col = "red")


str(preds.train <- predict(test, X, ind.row = ind.train, covar.row = covar[ind.train, ]))
mod <- ridge::linearRidge((y[ind.train] - rowMeans(preds.train)) ~ ., 
                          data = as.data.frame(preds.train), nPCs = 1)
summary(mod)

plot(apply(preds, 2, function(pred) cor(pred, y[ind.test])), 
     coefficients(mod)[-1])

preds3 <- preds2 + predict(mod, as.data.frame(preds))
plot(preds3, y[ind.test], pch = 20); abline(0, 1, col = "red")

cor(preds2, y[ind.test])
cor(preds3, y[ind.test])
