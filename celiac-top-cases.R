theone <- which.min(predict(gwas.train.gc))
table(G[ind.test, theone], y[ind.test])

counts.case <- big_counts(G, ind.row = which(y == 1))
counts.control <- big_counts(G, ind.row = which(y == 0))

counts.case[, theone]
counts.control[, theone]

props <- sapply(cols_along(G), function(j) {
  ccoj <- counts.control[, j]
  mut <- `if`(ccoj[1] < ccoj[3], 1, 3)
  `if`(counts.case[mut, j] > 20, counts.case[mut, j] / ccoj[mut], 0)
})
hist(props)
hist(props[props > 1])
ord.prop <- order(props, decreasing = TRUE)
head(ord.prop, 10)
lapply(ord.prop[1:10], function(j) table(G[, j], y))


# one geno
ind0 <- which(G[ind.test, ord.prop[1]] == 0)
mean(y[ind.test[ind0]])    ## 78.2%
# prs max
(M <- length(ind0))
indPRS <- order(prs[, which.max(aucs)], decreasing = TRUE)[1:M]
mean(y[ind.test[indPRS]])  ## 75.2%
# logit-simple (CMSA)
indLS <- order(preds2, decreasing = TRUE)[1:M]
mean(y[ind.test[indLS]])   ## 79.4 %
# logit-triple (CMSA)
indLT <- order(preds4, decreasing = TRUE)[1:M]
mean(y[ind.test[indLT]])   ## 82.4 %

plot(pROC::roc(y[ind.test], prs[, which.max(aucs)]))
plot(pROC::roc(y[ind.test], preds2), add = TRUE, col = "blue")
plot(pROC::roc(y[ind.test], preds4), add = TRUE, col = "red")

library(tidyverse)
library(plotROC)
scores_tidy <- tibble(
  d = y[ind.test], 
  PRS_max =  prs[, which.max(aucs)],
  LS = preds2,
  LT = preds4
) %>%
  gather(key = "Method", value = "Score", -d)
scores_tidy
  
bigstatsr:::MY_THEME(
  ggplot(scores_tidy, aes(d = d, m = Score, color = Method, linetype = Method)) +
    style_roc(xlab = "1 - Specificity", ylab = "Sensitivity")
) +
  geom_roc(n.cuts = 0, size = 2) +
  theme(legend.position = c(0.7, 0.3), legend.key.width = unit(4, "line")) +
  coord_equal()
  
  
