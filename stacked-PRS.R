library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
y <- celiac$fam$affection - 1L

U <- predict(readRDS("backingfiles/PCA.rds"))

set.seed(1)
ind.train <- sample(nrow(G), 5e3)
ind.test <- setdiff(rows_along(G), ind.train)

system.time(
  gwas <- big_univLogReg(G, y[ind.train], ind.train, covar.train = U[ind.train, ], 
                         ncores = nb_cores())
) # 140 sec


THR_P <- seq(0, max(-predict(gwas)))
THR_R2 <- c(0.01, 0.02, 0.05, 0.2, 0.5)

library(foreach)
test <- foreach(thr.r2 = THR_R2, .combine = "cbind") %do% {
  ind.keep <- snp_clumping(G, CHR, ind.row = ind.train,
                           S = abs(gwas$score), 
                           thr.r2 = thr.r2,
                           ncores = nb_cores())
  prs <- snp_PRS(G, betas.keep = gwas$estim[ind.keep],
                 ind.keep = ind.keep,
                 lpS.keep = -predict(gwas)[ind.keep],
                 thr.list = THR_P)
}
dim(test)
aucs <- apply(test[ind.test, ], 2, AUC, target = y[ind.test])
plot(aucs, pch = 20)
max(aucs)

plot(rep(THR_P, 5), aucs, pch = 20, col = rep(1:5, each = length(THR)), log = "x")

glm <- glmnet::glmnet(test[ind.train, ], y[ind.train], alpha = 1e-6, family = "binomial",
                      lower.limits = 0)
pred <- predict(glm, test[ind.test, ])
aucs.stacked <- apply(pred, 2, AUC, target = y[ind.test])
plot(aucs.stacked[-1], pch = 20)
max(aucs.stacked)



prs.train <- snp_PRS(G, betas.keep = gwas$estim[ind.keep],
                     ind.test = ind.train, 
                     ind.keep = ind.keep,
                     lpS.keep = -predict(gwas)[ind.keep],
                     thr.list = THR)
glm <- glmnet::glmnet(prs.train, y[ind.train], alpha = 0.1, family = "binomial",
                      lower.limits = 0)
aucs2 <- apply(predict(glm, prs), 2, AUC, target = y[ind.test])

glm <- glmnet::glmnet(prs.train, y[ind.train], alpha = 1, family = "binomial",
                      lower.limits = 0)
aucs3 <- apply(predict(glm, prs), 2, AUC, target = y[ind.test])

glm <- glmnet::glmnet(prs.train, y[ind.train], alpha = 0.001, family = "binomial",
                      lower.limits = 0)
aucs4 <- apply(predict(glm, prs), 2, AUC, target = y[ind.test])

# 0.05 / 0.2 / 0.5
max(aucs1)                             # 82.5 / 84.1 / 83.5
plot(aucs1, pch = 20)

max(aucs2)                             # 83.4 / 84.5 / 83.6
plot(aucs2[-1], pch = 20, col = "red")

max(aucs3)                             # 74.7 / 81.6 / 82.1
plot(aucs3[-1], pch = 20, col = "red")

max(aucs4)                             # 83.5 / 84.6 / 83.6
plot(aucs4[-1], pch = 20, col = "red")

# library(Matrix)
# rowSums(glm$beta)
# 
# library(ggplot2)
# qplot(1 / (1 + exp(-predict(glm, prs)[, 70])), fill = as.factor(y[ind.test]), 
#       alpha = I(0.4),  geom = "density") + 
#   theme_bigstatsr() +
#   scale_x_continuous(limits = c(0, 1))
