################################################################################

utils::globalVariables(
  c("par.causal", "par.dist", "par.model", "method", "num.simu",
    "true_set", "set", "pheno", "pred")
)

################################################################################

#' Alternative sensitivity & False Discovery Proportion
#'
#' Compute alternative sensitivity & False Discovery Proportion from
#' https://doi.org/10.1289/EHP172.
#'
#' @param x Data frame containing variables "true_set" and "set".
#' @param corr Sparse correlation matrix of the genotype matrix.
#'
#' @return `x` with 3 new variables: "AltSens", "AltFDP" and "nb.preds".
#' @import Matrix
#'
add_sens_FDP <- function(x, corr) {

  stopifnot(class(corr) == "dsCMatrix")
  corr <- methods::as(abs(corr), "dgTMatrix")

  alt <- purrr::map2(x$true_set, x$set, ~ {
    if (length(.y) == 0) return(c(0, 0))
    bool.x <- bool.y <- logical(nrow(corr))
    bool.x[.x] <- TRUE
    bool.y[.y] <- TRUE
    ind <- which_cond(corr@i, corr@j, bool.x, bool.y)
    AltSens <-     sum(tapply(corr@x[ind], corr@i[ind], max)) / length(.x)
    AltFDP  <- 1 - sum(tapply(corr@x[ind], corr@j[ind], max)) / length(.y)
    c(AltSens, AltFDP)
  }) %>%
    purrr::transpose()

  x %>%
    dplyr::mutate(AltSens = unlist(alt[[1]]),
                  AltFDP  = unlist(alt[[2]]),
                  nb.preds = purrr::map_int(set, length))
}

################################################################################

pAUC <- function(pred, target, p = 0.1) {
  val.min <- min(target)
  q <- stats::quantile(pred[target == val.min], probs = 1 - p)
  ind <- (target != val.min) | (pred > q)
  tryCatch(AUC(pred[ind], target[ind]) * p,
           error = function(e) p^2 / 2)  # all same prediction
}

add_AUC <- function(x, pAUC) {
  x %>%
    dplyr::mutate(AUC  = purrr::map2_dbl(pred, pheno, AUC),
                  pAUC = purrr::map2_dbl(pred, pheno, pAUC))
}

################################################################################

sub_list <- function(x, from, to) {
  for (i in seq_along(from)) {
    x <- sub(from[[i]], to[[i]], x)
  }
  x
}

rename_lvl <- function(x) {

  methods  <- c("PRS-all", "PRS-stringent", "PRS-max", "logit-simple", "logit-triple")
  methods2 <- c("C+T-all", "C+T-stringent", "C+T-max", "PLR",          "PLR3")

  dplyr::mutate(
    x,
    par.causal = factor(purrr::map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
    par.dist = stringr::str_to_title(par.dist),
    par.model = sub_list(par.model, c("simple", "fancy"), c("ADD", "COMP")),
    method = sub_list(method, methods, methods2)
  )
}

################################################################################

#' Read & format results
#'
#' @param files Files to read from.
#' @param corr Sparse correlation matrix of the genotype matrix.
#' @param ncores Number of cores to use. Default is `bigstatsr::nb_cores()`.
#'
#' @return A data frame of formated and aggregated results.
#' @export
#'
read_format_results <- function(files, corr, ncores = nb_cores()) {

  bigstatsr::big_apply(files, a.FUN = function(files, ind, corr) {

    files[ind] %>%
      purrr::map_dfr(~readRDS(.x)) %>%
      tibble::as_tibble() %>%
      rename_lvl() %>%
      add_AUC(pAUC) %>%
      add_sens_FDP(corr) %>%
      dplyr::select(-c(num.simu, true_set, set, pheno, pred))

  }, a.combine = bind_rows, ind = seq_along(files),
  ncores = ncores, block.size = 100, corr = corr)
}

################################################################################

#' Combine "method" and "r2"
#'
#' @param method Vector of names of methods.
#' @param r2 Vector of r2 values (possibly with missing values).
#'
#' @return Concatenation of `method` and `r2` using `"-"`.
#' @export
#'
#' @examples
#' c_method_r2(c("A", "B", "A"), c(0.5, NA, 0.2))
c_method_r2 <- function(method, r2) {
  ifelse(is.na(r2), method, paste0(method, "-", r2))
}

################################################################################

#' Bootstrap
#'
#' @param x Vector of results.
#' @param n Number of bootstrap replicate.
#' @param f Function to apply to samples of `x`. Default is `mean`.
#'
#' @return Standard deviation of bootstrap results.
#' @export
#'
boot <- function(x, n = getOption("nboot"), f = mean) {
  stats::sd(replicate(n, f(sample(x, replace = TRUE))))
}

################################################################################
