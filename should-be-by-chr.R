X <- matrix(rnorm(200), 20)
K <- solve(cor(X))
X_n <- scale(X) / sqrt(nrow(X) - 1)
all.equal(crossprod(X_n), cor(X))

svd <- svd(X_n, nu = 10, nv = 10)
mid <- tcrossprod(svd$u %*% diag(1 / svd$d^2))
K2 <- crossprod(X_n, mid) %*% X_n
all.equal(K, K2)


round(crossprod(svd$v), 4)
round(tcrossprod(svd$v), 4)

X <- matrix(rnorm(200), 10)
try(K <- solve(cor(X)))
X_n <- scale(X) / sqrt(nrow(X) - 1)
all.equal(crossprod(X_n), cor(X))

svd <- svd(X_n, nu = 9, nv = 9)
mid <- tcrossprod(sweep(svd$u, 2, svd$d[1:9]^2, '/'))
K2 <- crossprod(X_n, mid) %*% X_n

beta <- gwas.train$estim

# For celiac
cor.scaling <- function(X, ind.row, ind.col) {
  ms <- big_scale(center = TRUE, scale = TRUE)(X, ind.row, ind.col)
  ms$scale <- ms$scale * sqrt(length(ind.row) - 1)
  ms
}
system.time(
  K <- big_tcrossprodSelf(G, fun.scaling = cor.scaling, ind.row = ind.test)
) # 135 sec
mc <- attr(K, "center")
sc <- attr(K, "scale")

system.time(eigs <- eigen(K[], symmetric = TRUE))  ## 11 sec

str(eigs)
tail(eigs$values)

# ind <- head(seq_along(eigs$values), -1)
ind <- which(eigs$values > sqrt(ncol(K)))
mid <- tcrossprod(sweep(eigs$vectors[, ind], 2, eigs$values[ind], '/'))
y2 <- beta / sc
mid_beta <- mid %*% (big_prodVec(G, y2, ind.test) - drop(crossprod(y2, mc)))
pred <- big_prodVec(K, mid_beta)
AUC(pred, y[ind.test])
beta2 <- big_cprodVec(G, scale(mid_beta, scale = FALSE), ind.test) / sc

round(as.vector(by(cbind(beta, beta2), CHR, function(betas) cor(betas[, 1], betas[, 2]))), 3)

plot_density_AUC <- function(pred, true = y[ind.test]) {
  true_fct <- forcats::fct_recode(as.factor(true), Control = "0", Case = "1")
  data.frame(Score = pred, Status = true_fct) %>%
    ggplot() +
    geom_density(aes(Score, fill = Status), alpha = 0.4) +
    ggtitle(sprintf("AUC = %s%% | pAUC = %s", 
                    round(100 * AUC(pred, true), 2),
                    round(pAUC(pred, true), 4))) + 
    theme_bigstatsr() +
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
}
(p1 <- plot_density_AUC(pred))


system.time(
  K <- big_tcrossprodSelf(G, fun.scaling = cor.scaling, ind.test)
) # 135 sec
system.time(eigs <- eigen(K[]))  ## 11 sec

# ind <- head(seq_along(eigs$values), -1)
ind <- which(eigs$values > sqrt(ncol(K)))
mid <- tcrossprod(sweep(eigs$vectors[, ind], 2, eigs$values[ind], '/'))

# Chromosome 6 only
test <- lapply(1:22, function(chr) {
  print(chr)
  ind.chr <- which(CHR == chr)
  
  mc <- attr(K, "center")[ind.chr]
  sc <- attr(K, "scale")[ind.chr]
  
  y2 <- beta[ind.chr] / sc
  mid_beta <- mid %*% scale(big_prodVec(G, y2, ind.test, ind.chr), scale = FALSE)
  pred <- big_prodVec(K, mid_beta)
  beta2 <- big_cprodVec(G, scale(mid_beta, scale = FALSE), ind.test, ind.chr) / sc
  list(pred = pred, beta = beta2)
})

sapply(test, function(res) AUC(res$pred, y[ind.test]))
AUC(rowSums(sapply(test, function(res) res$pred)), y[ind.test])

plot(beta[CHR == chr], beta2, pch = 20)

### K chr by chr does not work
# [1] 0.5285619 0.5267598 0.5337792 0.4877519 0.5086177 0.7339950 0.5205857 0.5187503
# [9] 0.5228855 0.4968187 0.5215254 0.5029465 0.5071346 0.5088806 0.5159409 0.5296606
# [17] 0.5112856 0.5183893 0.5041300 0.4967401 0.4970916 0.4998657

### K for all
# [1] 0.5377430 0.5522180 0.5114261 0.5328115 0.5059283 0.7172683 0.5439647 0.5295968
# [9] 0.5006401 0.5057092 0.4960067 0.5287786 0.4867427 0.5094221 0.4941526 0.5324386
# [17] 0.5498858 0.5146355 0.4965682 0.4946375 0.5139183 0.5004015