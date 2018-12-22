################################################################################


#' PRS (Clumping + Thresholding)
#'
#' @inheritParams bigsnpr::snp_clumping
#' @param pheno.all All phenotypes.
#' @param covar.all All covariables.
#' @param ind.train Indices corresponding to the training set.
#' @param ind.test Indices corresponding to the test set.
#' @param ncores Number of cores to use. Default is `bigstatsr::nb_cores()`.
#' @param family Default is `"logistic"`. Can also use `"linear"`.
#'
#' @return A data frame with 5 variables:
#'   - "method": the name of the method,
#'   - "pred": the predictions for the test set,
#'   - "thr.r2": the threshold used for clumping,
#'   - "set": the set of non-zero coefficients,
#'   - "timing": the execution time of the GWAS step only.
#' @export
#' @import bigsnpr
#'
PRS <- function(G, infos.chr, infos.pos,
                pheno.all, covar.all, ind.train, ind.test,
                ncores = nb_cores(), family = "logistic") {

  fun_mod <- `if`(family == "logistic", bigstatsr::big_univLogReg,
                  bigstatsr::big_univLinReg)
  fun_msr <- `if`(family == "logistic", bigstatsr::AUC, stats::cor)

  timing <- system.time({

    # GWAS
    gwas.train <- fun_mod(
      G, pheno.all[ind.train], ind.train = ind.train,
      covar.train = covar.all[ind.train, , drop = FALSE],
      ncores = ncores
    )
    gwas.train.gc <- snp_gc(gwas.train)

  })[3]

  lapply(c(0.05, 0.2, 0.8), function(thr.r2) {

    # Clumping
    ind.keep <- snp_clumping(G, infos.chr = infos.chr,
                             ind.row = ind.test,
                             thr.r2 = thr.r2,
                             S = abs(gwas.train.gc$score),
                             size = 500,
                             is.size.in.bp = TRUE,
                             infos.pos = infos.pos,
                             ncores = ncores)

    # PRS
    thrs <- c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100)))
    lpS <- -stats::predict(gwas.train.gc)
    prs <- snp_PRS(G, betas.keep = gwas.train.gc$estim[ind.keep],
                   ind.test = ind.test,
                   ind.keep = ind.keep,
                   lpS.keep = lpS[ind.keep],
                   thr.list = thrs)

    ind.best <- which.max(apply(prs, 2, fun_msr, target = pheno.all[ind.test]))

    methods <- c("PRS-all", "PRS-stringent", "PRS-max")
    indices <- c(1:2, ind.best)

    lapply(1:3, function(i) {
      k <- indices[i]
      tibble(
        method = methods[i],
        pred = list(prs[, k]),
        thr.r2 = thr.r2,
        set = list(intersect(ind.keep, which(lpS > thrs[k]))),
        timing = timing
      )
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}

################################################################################

#' Sparse Logistic Regression
#'
#' @inheritParams PRS
#' @param method Name of the method.
#' @param alphas See [bigstatsr::big_spLogReg].
#'
#' @export
#' @return A data frame with 5 variables:
#'   - "method": the name of the method,
#'   - "pred": the predictions for the test set,
#'   - "timing": the execution time,
#'   - "alpha": the alpha that maximized prediction for the validation sets,
#'   - "set": the set of non-zero coefficients.
#'
logit.CMSA <- function(G, pheno.all, covar.all, ind.train, ind.test, method,
                       alphas = c(1, 0.5, 0.05, 0.001),
                       ncores = nb_cores(), family = "logistic") {

  timing <- system.time({

    cmsa.logit <- `if`(
      family == "logistic",
      bigstatsr::big_spLogReg,
      bigstatsr::big_spLinReg
    )(
      X = G, y01.train = pheno.all[ind.train],
      ind.train = ind.train,
      covar.train = covar.all[ind.train, , drop = FALSE],
      alphas = alphas, ncores = ncores
    )

    preds <- stats::predict(cmsa.logit, X = G, ind.row = ind.test,
                            covar.row = covar.all[ind.test, , drop = FALSE])
  })[3]

  tibble(
    method   = method,
    pred     = list(preds),
    timing   = timing,
    alpha = attr(cmsa.logit, "alpha"),
    set = list(which(rowSums(sapply(cmsa.logit, function(x) x$beta.X) != 0) > 0))
  )
}

################################################################################

#' T-Trees
#'
#' @param TTree Path of TTree executable.
#' @param file.base Path of file containing the genotype matrix.
#' @inheritParams PRS
#' @param n.trees Number of trees.
#'
#' @return A data frame with 5 variables:
#'   - "method": the name of the method,
#'   - "pred": the predictions for the test set,
#'   - "timing": the execution time,
#'   - "set": the set of non-zero coefficients.
#' @export
#' @importFrom utils read.table
#'
ttrees <- function(TTree, file.base, pheno.all, ind.train, ind.test,
                   n.trees = 100) {

  stopifnot(file.exists(TTree))

  file0.jdb  <- paste0(file.base, ".jdb")
  file0.bloc <- paste0(file.base, ".bloc")

  # Write jdb file with new pheno
  tmpfile <- tempfile()
  write(c(rep("", 5), pheno.all), ncolumns = 1, tmpfile)
  # https://stackoverflow.com/a/7846550/6103040
  system(sprintf("awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' %s %s > %s",
                 tmpfile, file0.jdb, file.jdb <- paste0(tmpfile, ".jdb")))
  # Write indices of learning and validation sets
  file.learn <- paste0(tmpfile, "_learn.txt")
  file.val   <- paste0(tmpfile, "_val.txt")
  cat(ind.train - 1, file = file.learn, sep = "\t")
  cat(ind.test - 1,  file = file.val,   sep = "\t")

  timing <- system.time(
    system(glue::glue(
      "{TTree}",
      " -j {file.jdb}",
      " -m 3",
      " -b {file0.bloc}",
      " -l {file.learn}",
      " -v {file.val}",
      " -t {n.trees} -k 1000 -c 5 -n 2000",
      " -x -s"
    ))
  )[3]

  file.roc <- sprintf("%s_k1000_m3_t%d_ic5_nmin2000_0000.roc",
                      file.jdb, n.trees)
  file.vim <- sub("\\.roc$", ".vim", file.roc)
  preds <- read.table(file.roc, header = FALSE)

  tibble(
    method   = "T-Trees",
    pred     = list(preds[match(ind.test - 1, preds[[1]]), 2]),
    timing   = timing,
    set = list(which(read.table(file.vim)$V2 != 0))
  )
}

################################################################################
