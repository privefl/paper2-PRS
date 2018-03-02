N <- 1000
x1 <- rnorm(N)
x2 <- rnorm(N, mean = 1)
x <- c(x1, x2)
y <- rep(0:1, each = N)
plot(print(pROC::roc(y ~ x)))
bigstatsr::AUC(x, y)
plot(print(pROC::roc(y[ind] ~ x[ind])))
bigstatsr::AUC(x[x > q], y[x > q])
plot(print(pROC::roc(y[x > q] ~ x[x > q])))


test <- rocdemo.sca(y, x)
ROC::AUC(test)
plot(test)
sapply(0:10 / 10, function(p) ROC::pAUC(test, t0 = p))
plot(test)

# Partial AUC based on AUC
pAUC <- function(pred, target, p = 0.1) {
  q <- quantile(pred[target == 0], probs = 1 - p)
  ind <- (target > 0.5) | (pred > q)
  bigstatsr::AUC(pred[ind], target[ind]) * p
}
sapply(1:9 / 10, function(p) pAUC(x, y, p))

plot(pROC::roc(y ~ x))
abline(v = 0.9, col = "red")
abline(h = mean(x2 > quantile(x1, 0.9)), col = "blue")
