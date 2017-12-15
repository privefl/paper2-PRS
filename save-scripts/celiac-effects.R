# Effects GWAS vs logistic regression
qplot(cmsa.logit[cols_along(G)], gwas.train.gc$estim,
      color = predict(gwas.train.gc, log10 = FALSE) < 5e-8) %>%
  bigstatsr:::MY_THEME() +
  labs(color = "Bonf") + 
  geom_abline(slope = 1) + 
  geom_smooth(method = "lm", color = "green")


preds2


# obj.svd3 <- obj.svd
# obj.svd3$u[ind.train, ] <- NA
# plot(obj.svd3, type = "scores")
# 
# 
# preds.train <- predict(cmsa.logit, X = G, ind.row = ind.train, 
#                        covar.row = obj.svd$u[ind.train, ])
# myglm <- glm(y[ind.train] ~ preds.train, family = "binomial")
# summary(myglm)
# preds2.prob <- predict(myglm, newdata = data.frame(preds.train = preds2),
#                        type = "response")
# err <- rep(NA, nrow(G))
# err[ind.test] <- (y[ind.test] != (preds2.prob > 0.5)) # (y[ind.test] - preds2.prob)^2
# 
# plot(obj.svd3, type = "scores") + 
#   aes(color = as.factor(err), alpha = I(rep(0.6, nrow(G))))

boot1000 <- function(x) {
  sd(replicate(1000, mean(sample(x, replace = TRUE))))
}
round4 <- function(x) round(x, 4)
data.frame(pop = celiac$fam$family.ID[ind.test],
           pred = preds2.prob,
           error = err[ind.test]) %>%
  filter(pred > quantile(pred, 0.8)) %>%
  group_by(pop) %>%
  summarise_at("error", funs(mean, boot1000)) %>%
  mutate_if(is.numeric, funs(round4))


head(ord.eff <- order(abs(cmsa.logit[cols_along(G)]), decreasing = TRUE), n = 10)
lapply(ord.eff[1:10], function(j) table(G[, j], y))
# Manhattan plots with greatest effects
cowplot::plot_grid(
  snp_manhattan(gwas.train.gc, CHR, POS, labels = labels, 
                ind.highlight = head(ord.eff, 50)),
  snp_manhattan(gwas.train.gc, CHR, POS, labels = labels, 
                ind.highlight = head(ord.eff, 50)) +
    coord_cartesian(ylim = c(0, 25)) + 
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 3),
  align = "hv", ncol = 1, labels = LETTERS[1:2], label_size = 25, scale = 0.95
)
ggsave("colored-manhattan.png", scale = 1/90, width = 1270, height = 1040)

# Is the first effects correlated? Osef
round(cor(G[ind.train, head(ord.eff, 10)])^2, 2)
cmsa.logit[head(ord.eff, 10)]


data.frame(pop = celiac$fam$family.ID[ind.test],
           pred = preds2,
           true = y[ind.test],
           geno = G[ind.test, ord.eff[3]]) %>%
  myggplot() + 
  geom_density(aes(pred, fill = as.factor(true)), alpha = 0.3) +
  facet_grid(pop ~ factor(geno, levels = 2:0)) +
  labs(fill = "Pheno")


