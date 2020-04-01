################################################################################

scaled_prod <- function(X, beta, ind.col) {
  ms <- big_scale()(X, rows_along(X), ind.col)
  big_prodVec(X, beta, ind.col = ind.col, center = ms$center, scale = ms$scale)
}

#' Simulate phenotypes
#'
#' @param G Matrix of genotypes.
#' @param h2 Heritability.
#' @param M Number of causal variants. Should be multiple of 6.
#' @param ind.possible Indices of possible causal variants.
#' @param effects.dist Distribution of effects.
#'   Either `"gaussian"` (the default) or `"laplace"`.
#' @param model Model for simulation.
#'   Either `"simple"` (the default) with linear effets only,
#'   or `"fancy"` with linear, domination and interaction-type effetcs.
#' @param K Prevalence. Default is `0.3`.
#'
#' @return A list with 3 elements:
#'   - "pheno": vector of binary phenotypes,
#'   - "set": indices of causal variants,
#'   - "effects": Size of effects corresponding to `set`.
#' @export
#' @importFrom stats qnorm rnorm sd var
#'
#' @import bigstatsr
#'
get_pheno <- function(
  G,                                        ## matrix of genotypes
  h2,                                       ## heritability
  M,                                        ## nbs of causal variants
  ind.possible = bigstatsr::cols_along(G),  ## indices of possible causal variants
  effects.dist = c("gaussian", "laplace"),  ## distribution of effects
  model = c("simple", "fancy"),             ## model for simulation
  K = 0.3                                   ## prevalence
) {

  effects.dist  <- match.arg(effects.dist)
  model <- match.arg(model)

  set <- sample(ind.possible, size = M)
  effects <- `if`(effects.dist == "gaussian",
                  rnorm(M, sd = sqrt(h2 / M)),
                  rmutil::rlaplace(M, s = sqrt(h2 / (2*M))))

  if (model == "simple") {
    # only linear
    y.simu <- scaled_prod(G, effects, set)
  } else {
    sets <- split(1:M, sample(rep_len(1:3, M)))
    # linear
    ind1 <- sets[[1]]
    y.simu <- scaled_prod(G, effects[ind1], set[ind1])
    # recessive / dominant
    ind2 <- sets[[2]]
    y.simu <- y.simu + scaled_prod(G$copy(code = (G$code256 > 0.5)),
                                   effects[ind2], set[ind2])
    # interactions
    ind3 <- matrix(sets[[3]], ncol = 2)
    y.simu <- y.simu + scale(G[, set[ind3[, 1]]] * G[, set[ind3[, 2]]]) %*%
      effects[ind3[, 1]] * sqrt(2)
  }

  y.simu <- y.simu / sd(y.simu) * sqrt(h2)
  stopifnot(all.equal(drop(var(y.simu)), h2))
  y.simu <- y.simu + rnorm(nrow(G), sd = sqrt(1 - h2))
  list(pheno = as.integer(y.simu > qnorm(1 - K)), set = set, effects = effects)
}

################################################################################
