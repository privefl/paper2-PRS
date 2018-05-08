library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
hapmap3 <- snp_attach("../2017_project_5/data/hapmap3_qc.rds")

# imputation
rds <- "../2017_project_5/data/hapmap3_qc-infos-impute.rds"
tmp <- big_attach(rds)

hasNA <- is.na(tmp[1, ])
ind <- match(hapmap3$map$marker.ID, celiac$map$marker.ID)
obj.svd0 <- snp_autoSVD(
  G = hapmap3$genotypes$copy(bigsnpr:::CODE_IMPUTE_PRED),
  infos.chr = hapmap3$map$chromosome,
  infos.pos = hapmap3$map$physical.pos,
  ind.col = which(!is.na(ind) & !hasNA), 
  k = 20, 
  ncores = nb_cores()
)

plot(obj.svd0) + 
  ggplot2::scale_y_log10()
plot(obj.svd0, type = "loadings", loadings = 1:6, coeff = 0.5)
plot(obj.svd0, type = "scores", scores = 9:10)

test <- varimax(obj.svd0$v)

G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
obj.svd <- readRDS("backingfiles/PCA.rds")
u <- obj.svd$u

u0 <- predict(obj.svd0, G, ind.col = ind[attr(obj.svd0, "subset")])
u0_rot <- u0 %*% test$rotmat

files <- list.files("backingfiles/others/", "^cluster",
                    full.names = TRUE)
pop <- snp_getSampleInfos(celiac, files)[[1]]
pop[pop == 4] <- 3

plot(u0, col = pop, pch = 20)
plot(u0_rot, col = pop, pch = 20)
plot(u0[, 3:4], col = pop, pch = 20)
plot(u0_rot[, 3:4], col = pop, pch = 20)


# PRS
y <- celiac$fam$affection - 1
ind.train <- sort(sample(nrow(G), 10e3))
ind.test <- setdiff(rows_along(G), ind.train)
system.time(
  prs0 <- big_spLogReg(G, y01.train = y[ind.train], 
                       ind.train = ind.train,
                       covar.train = u[ind.train, ], 
                       ncores = nb_cores())
) # 170


system.time(
  prs5 <- big_spLogReg(G, y01.train = y[ind.train], 
                       ind.train = ind.train,
                       covar.train = cbind(u0, u0_rot)[ind.train, ], 
                       ncores = nb_cores())
) # 182

system.time(
  prs05 <- big_spLogReg(G, y01.train = y[ind.train], 
                        ind.train = ind.train,
                        covar.train = cbind(u, u_rot5)[ind.train, ], 
                        ncores = nb_cores())
) # 182

pred0 <- rowMeans(predict(prs0, G, ind.test, covar.row = u[ind.test, ], proba = FALSE))
AUC(pred0, y[ind.test])
pred5 <- rowMeans(predict(prs5, G, ind.test, covar.row = cbind(u0, u0_rot)[ind.test, ], proba = FALSE))
AUC(pred5, y[ind.test])
pred05 <- rowMeans(predict(prs05, G, ind.test, covar.row = cbind(u, u_rot5)[ind.test, ], proba = FALSE))
AUC(pred05, y[ind.test])

# str(prs0)
sapply(prs0, function(x) x$beta.covar)

# str(prs5)
sapply(prs5, function(x) x$beta.covar)

preds <- data.frame(pred0 = pred0, 
                    pred5 = pred5,
                    # pred05 = pred05,
                    true  = y[ind.test],
                    pop   = pop[ind.test])
library(ggplot2)
qplot(pred0, geom = "density", data = preds,
      fill = as.factor(pop),
      alpha = I(0.4)) + 
  facet_grid(~ true)

qplot(pred5, geom = "density", data = preds,
      fill = as.factor(pop),
      alpha = I(0.4)) + 
  facet_grid(~ true)

qplot(pred05, geom = "density", data = preds,
      fill = as.factor(pop),
      alpha = I(0.4)) + 
  facet_grid(~ true)

library(dplyr)
preds %>%
  group_by(true, pop) %>%
  summarise_all(c(mean, sd))