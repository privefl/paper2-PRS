plink2 <- bigsnpr::download_plink2("tmp-data")

ids <- bigsnpr::snp_attach("data/ukbb.rds")$fam[1:2]

if (!file.exists("data/ukbb.pgen")) {
  bigsnpr:::write.table2(ids, tmp <- tempfile())
  system(glue::glue("{plink2} --bfile data/ukbb --keep {tmp}",
                    " --make-pgen vzs --out data/ukbb",
                    " --threads 16 --memory 200000"))
}
file.size("data/ukbb.pgen") / 1024^3  # 22 GB

# Make sure individuals are in the sale order
library(bigreadr)
all.equal(fread2("data/ukbb.psam", select = 1:2), ids, check.attributes = FALSE)

library(dplyr)
df0 <- readRDS("data/pheno_covar.rds")
covar <- select(df0, sex, age, PC1:PC16)
pheno <- select(df0, height, BMI, asthma, high_cholesterol)
set <- df0$set
split <- case_when(set == 5 ~ "test", set == 1 ~ "val", TRUE ~ "train")
table(split)
#  test  train    val
# 67495 202485  67495

phenotype.file <- bigreadr::fwrite2(cbind.data.frame(
  FID = ids[[1]],
  IID = ids[[2]],
  split,
  split_refit = case_when(set == 5 ~ "val", TRUE ~ "train"),
  pheno,
  covar
), "data/pheno_snpnet.txt")

bigassertr::assert_dir("snpnet2")
bigassertr::assert_dir("log")

NCORES <- 16
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "2-00:00", c = NCORES, mem = "500g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_iwalk(pheno, ~ {

  res_file <- paste0("snpnet2/", .y, ".rds")
  if (file.exists(res_file)) return(NULL)

  unlink(res_dir <- paste0("tmp-data/", .y, 2), recursive = TRUE)

  configs <- list(plink2.path = plink2, nCores = NCORES, niter = 100,
                  use.glmnetPlus = TRUE, results.dir = res_dir,
                  early.stopping = TRUE, stopping.lag = 2)

  ## Penalized reg
  library(snpnet)
  time <- system.time(
    fit_snpnet <- snpnet(
      genotype.pfile = "data/ukbb",
      phenotype.file = phenotype.file,
      phenotype = .y,
      covariates = names(covar),
      # alpha = 1,
      family = if (all(.x %in% c(0:1, NA))) "binomial" else "gaussian",
      split.col = "split",
      # mem = 200e3,
      configs = configs
    )
  )[3]

  saveRDS(list(mod = fit_snpnet, time = time), res_file)
})

# Refitting
furrr::future_iwalk(pheno, ~ {

  res_file <- paste0("snpnet2/", .y, "2.rds")
  if (file.exists(res_file)) return(NULL)

  fit_snpnet_train <- readRDS(paste0("snpnet2/", .y, ".rds"))$mod
  opt_idx <- which.max(fit_snpnet_train$metric.val)

  unlink(res_dir <- paste0("tmp-data/", .y, 3), recursive = TRUE)

  configs <- list(plink2.path = plink2, nCores = NCORES, niter = 100,
                  use.glmnetPlus = TRUE, results.dir = res_dir,
                  early.stopping = FALSE)

  ## Penalized reg
  library(snpnet)
  time <- system.time(
    fit_snpnet <- snpnet(
      genotype.pfile = "data/ukbb",
      phenotype.file = phenotype.file,
      phenotype = .y,
      covariates = names(covar),
      lambda = fit_snpnet_train$full.lams[1:opt_idx],
      # alpha = 1,
      family = if (all(.x %in% c(0:1, NA))) "binomial" else "gaussian",
      split.col = "split_refit",
      # mem = 200e3,
      configs = configs
    )
  )[3]

  saveRDS(list(mod = fit_snpnet, time = time), res_file)
})
