
library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
obj.svd <- readRDS("backingfiles/PCA.rds")
plot(obj.svd) + 
  ggplot2::scale_y_log10(breaks = c(500, 1000, 2000))

files <- list.files("backingfiles/others/", "^cluster",
                    full.names = TRUE)
pop <- snp_getSampleInfos(celiac, files)[[1]]
pop[pop == 4] <- 3
library(ggplot2)
plot(obj.svd, type = "scores", scores = 3:4) +
  aes(color = as.factor(pop))


test <- varimax(tmp$u)
plot(test$loadings)
plot(test$loadings[, 3:4])


tmp2 <- prcomp(iris[1:4], center = TRUE, scale. = TRUE)
plot(tmp2$x, pch = 20, col = iris$Species)
rot_var <- varimax(tmp2$rotation[, 1:2])
rot_var
plot(rot_var$loadings)
plot(tmp2$x[, 1:2] %*% rot_var$rotmat, pch = 20, col = iris$Species)

# test on PCs of Celiac
test <- varimax(tmp$v)
test2 <- varimax(tmp$v[, 1:3])
test2$rotmat
test$rotmat[1:3, 1:3]
str(test)     
test$loadings
plot(test$rotmat)
test$rotmat
plot(tmp$u, pch = 20)
u_rot <- tmp$u %*% test$rotmat
u_rot3 <- tmp$u[, 1:3] %*% test2$rotmat
plot(u_rot, col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot3, col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot[, 3:4], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)
u <- obj.svd$u
plot(u[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)


test3 <- varimax(obj.svd$u, normalize = FALSE)
u_rot2 <- u %*% test3$rotmat
plot(u, col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot2, col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot2[, 3:4], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot2[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)


## Test promax
test4 <- promax(u)
all.equal(u_rot4 <- unclass(test4$loadings), u %*% test4$rotmat)
plot(u, col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot4, col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot4[, 3:4], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot4[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)

## Test promax on loadings
test5 <- promax(obj.svd$v)
# test5 <- varimax(obj.svd$v, normalize = TRUE, eps = 1e-10)
# test5$rotmat[1:6, 1:6]
# test6 <- promax(obj.svd$v[, 1:6])
# test6$rotmat
# 
# plot(test5$loadings[, 1], obj.svd$v[, 1], pch = 20)
# plot(test5$loadings[, 9], obj.svd$v[, 9], pch = 20)
# 
# bigstatsr:::predict.big_SVD

u_rot5 <- u %*% test5$rotmat
plot(u_rot5[, 1:2], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, 3:4], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, 7:8], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, 9:10], col = scales::alpha(pop, 0.2), pch = 20)
plot(u[, 9:10], col = scales::alpha(pop, 0.2), pch = 20)

u_rot6 <- u[, 1:6] %*% test6$rotmat
plot(u_rot5[, 1:6], u_rot6)

plot(svd(u_rot5)$u[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)
gwas0 <- snp_pcadapt(G, u)
snp_manhattan(gwas0, CHR, POS, npoints = 20e3)
gwas5 <- snp_pcadapt(G, u_rot5)
snp_manhattan(gwas5, CHR, POS, npoints = 20e3)
all.equal(gwas0$score, gwas5$score)
plot(-predict(gwas0), -predict(gwas5), pch = 20)
# SAME FOR PCADAPT

tmp_u <- varimax(u)
u_rot7 <- tmp_u$loadings
plot(u_rot7[, 1:2], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot7[, 3:4], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot7[, 5:6], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot7[, 7:8], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot7[, 9:10], col = scales::alpha(pop, 0.2), pch = 20)



# GWAS
y <- celiac$fam$affection - 1
ind6 <- which(CHR %in% 6)
gwas0 <- big_univLogReg(G, y, ind.col = ind6, covar.train = u,
                        ncores = nb_cores())
snp_qq(gwas0)
gwas5 <- big_univLogReg(G, y, ind.col = ind6, covar.train = u_rot5,
                        ncores = nb_cores())
snp_qq(gwas5)
all.equal(gwas0$score, gwas5$score)

# PRS
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
                       covar.train = u_rot5[ind.train, ], 
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
pred5 <- rowMeans(predict(prs5, G, ind.test, covar.row = u[ind.test, ], proba = FALSE))
AUC(pred5, y[ind.test])
pred05 <- rowMeans(predict(prs05, G, ind.test, covar.row = cbind(u, u_rot5)[ind.test, ], proba = FALSE))
AUC(pred05, y[ind.test])
preds <- data.frame(pred0 = pred0, 
                    pred5 = pred5,
                    pred05 = pred05,
                    true  = y[ind.test],
                    pop   = pop[ind.test])
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

# str(prs0)
sapply(prs0, function(x) x$beta.covar)

# str(prs5)
sapply(prs5, function(x) x$beta.covar)
sapply(prs05, function(x) x$beta.covar)

plot(u_rot5[, 3], u_rot5[, 4], col = scales::alpha(pop, 0.4), pch = 20)

predict(big_univLinReg(big_copy(cbind(u, u_rot5)), y), log10 = FALSE)

plot(u[, c(4, 6)], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, c(4, 6)], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, 3:4], col = scales::alpha(pop, 0.4), pch = 20)
plot(u_rot5[, 2:3], col = scales::alpha(pop, 0.4), pch = 20)



sex <- celiac$fam$sex
predict(big_univLinReg(big_copy(cbind(u, u_rot5)), sex), log10 = FALSE)
plot(u_rot5[, c(2, 10)], col = scales::alpha(sex, 0.4), pch = 20)
