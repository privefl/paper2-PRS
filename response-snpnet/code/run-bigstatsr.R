library(bigsnpr)
library(dplyr)

ukbb <- snp_attach("data/ukbb.rds")
G <- ukbb$genotypes
file.size(G$bk) / 1024^3  # 158 GB

df0 <- readRDS("data/pheno_covar.rds")
covar <- select(df0, sex, age, PC1:PC16)
pheno <- select(df0, height, BMI, asthma, high_cholesterol)
set <- df0$set

bigassertr::assert_dir("bigstatsr")
bigassertr::assert_dir("log")

library(future.batchtools)
NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "128g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_imap(pheno, ~ {

  res_file <- paste0("bigstatsr/", .y, ".rds")
  if (file.exists(res_file)) return(NULL)

  y <- .x
  ind.train <- which(!is.na(y) & complete.cases(covar) & set != 5)

  FUN <- if (all(y[ind.train] %in% 0:1)) big_spLogReg else big_spLinReg

  ## Penalized reg
  time <- system.time(
    mod <- FUN(G, y[ind.train], ind.train = ind.train,
               covar.train = as.matrix(covar[ind.train, ]),
               pf.covar = rep(0, ncol(covar)),
               alphas = 1,         # lasso
               power_scale = 0,    # not scaled
               ind.sets = set[ind.train], K = 4,
               n.abort = 2, nlam.min = 30,
               ncores = NCORES)
  )[3]

  saveRDS(list(mod = mod, time = time), res_file)
})
