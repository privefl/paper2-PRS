library(bigsnpr)
ukbb <- snp_attach("data/ukbb.rds")
G <- ukbb$genotypes

df0 <- readRDS("data/pheno_covar.rds")
covar <- dplyr::select(df0, sex, age, PC1:PC16)
ind.test <- which(df0$set == 5)

# perf <- function(pred, target) {
#   is_not_na <- which(!is.na(target))
#   AUC(pred[is_not_na], target[is_not_na])
# }

perf <- function(pred, target) {
  my_pcor <- function(x, y, z) {
    df <- cbind.data.frame(z, x = x, y = y)
    cor(lm(x ~ . - y, data = df)$residuals,
        lm(y ~ . - x, data = df)$residuals,
        use = "pairwise.complete.obs")
  }
  my_pcor(pred, target, covar[ind.test, ])
}

# name <- "asthma"
# name <- "high_cholesterol"

for (name in c("height", "BMI", "asthma", "high_cholesterol")) {

  res <- readRDS(paste0("bigstatsr/", name, ".rds"))
  y <- df0[[name]]

  beta_one <- head(res$mod[[1]][[1]]$beta, -ncol(covar))
  ind <- which(beta_one != 0)
  pred_one <- big_prodVec(G, beta_one[ind], ind.row = ind.test,
                          ind.col = ind, ncores = 15)

  pred_all <- predict(res$mod, G, ind.row = ind.test,
                      covar.row = matrix(0, length(ind.test), ncol(covar)),
                      proba = FALSE, ncores = 15)

  res2 <- readRDS(paste0("snpnet2/", name, ".rds"))
  id <- which.max(res2$mod$metric.val)
  beta_snpnet <- tail(res2$mod$beta[[id]], -ncol(covar))
  ind <- match(names(beta_snpnet), with(ukbb$map, paste0(marker.ID, "_", allele1)))
  pred_snpnet <- big_prodVec(G, beta_snpnet, ind.row = ind.test,
                             ind.col = ind, ncores = 15)

  res3 <- readRDS(paste0("snpnet2/", name, "2.rds"))
  id <- length(res3$mod$metric.train)
  beta_snpnet2 <- tail(res3$mod$beta[[id]], -ncol(covar))
  ind <- match(names(beta_snpnet2), with(ukbb$map, paste0(marker.ID, "_", allele1)))
  pred_snpnet2 <- big_prodVec(G, beta_snpnet2, ind.row = ind.test,
                              ind.col = ind, ncores = 15)

  print(list(
    trait = name,
    time = round(c(
      all_bigstatsr = res$time,
      one_snpnet = res2$time,
      all_snpnet = res3$time
    ) / 60),
    perf = round(c(
      one_bigstatsr = perf(pred_one, y[ind.test]),
      all_bigstatsr = perf(pred_all, y[ind.test]),
      one_snpnet = perf(pred_snpnet, y[ind.test]),
      all_snpnet = perf(pred_snpnet2, y[ind.test])
    ), 4)
  ))

}

library(ggplot2)
qplot(pred_one, pred_snpnet) +
  bigstatsr::theme_bigstatsr() +
  geom_abline(color = "chartreuse") +
  coord_equal() +
  labs(x = "Genetic scores from snpnet",
       y = "Genetic scores from bigstatsr")

qplot(pred_all, pred_snpnet2) +
  bigstatsr::theme_bigstatsr() +
  geom_abline(color = "chartreuse") +
  coord_equal() +
  labs(x = "Genetic scores from snpnet",
       y = "Genetic scores from bigstatsr")
