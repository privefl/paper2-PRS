
triple_biglasso <- function(X, y.train, ind.train, covar.train,
                            alpha = 1,
                            lambda.min = `if`(nrow(X) > ncol(X), .001, .01),
                            nlambda = 100, lambda.log.scale = TRUE,
                            lambda, eps = 1e-7, max.iter = 1000,
                            dfmax = 20e3,
                            penalty.factor = NULL,
                            warn = TRUE,
                            verbose = FALSE) {
  
  lambda.min <- max(lambda.min, 1e-6)
  
  n <- length(ind.train) ## subset of X. idx: indices of rows.
  if (is.null(covar.train)) covar.train <- matrix(0, n, 0)
  assert_lengths(y.train, ind.train, rows_along(covar.train))
  
  p <- 3 * ncol(X) + ncol(covar.train)
  if (is.null(penalty.factor)) penalty.factor <- rep(1, p)
  if (p != length(penalty.factor))
    stop("'penalty.factor' has an incompatible length.")
  
  if (alpha == 1) {
    penalty <- "lasso"
  } else if (alpha < 1 && alpha > 1e-6) {
    penalty <- "enet"
  } else {
    stop("alpha must be between 1e-6 and 1 for elastic net penalty.")
  }
  
  if (nlambda < 2) stop("nlambda must be at least 2")
  
  if (any(is.na(y.train)))
    stop(paste("Missing data (NA's) detected. Take actions",
               "(e.g., removing cases, removing features, imputation)",
               "to eliminate missing data before fitting the model."))
  
  if (class(y.train) != "numeric")
    tryCatch(y.train <- as.numeric(y.train), error = function(e)
      stop("y.train must numeric or able to be coerced to numeric"))
  
  yy <- transform_levels(y.train)
  
  if (missing(lambda)) {
    user.lambda <- FALSE
    lambda <- rep(0, nlambda);
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  ## fit model
  if (verbose) printf("\nStart biglasso: %s\n", format(Sys.time()))
  
  Rcpp::sourceCpp('biglassoLogTriple.cpp')
  res <- triple_cdfit_binomial_hsr(X, yy, ind.train, covar.train,
                                   lambda, nlambda, lambda.log.scale,
                                   lambda.min, alpha,
                                   user.lambda | any(penalty.factor == 0),
                                   eps, max.iter, penalty.factor,
                                   dfmax, warn, verbose)
  
  
  a <- res[[1]]
  b <- Matrix(res[[2]], sparse = TRUE)
  center <- res[[3]]
  scale <- res[[4]]
  lambda <- res[[5]]
  loss <- res[[6]]
  iter <- res[[7]]
  rejections <- res[[8]]
  col.idx <- res[[9]]
  
  if (verbose) printf("\nEnd biglasso: %s\n", format(Sys.time()))
  
  # p.keep <- length(col.idx)
  col.idx <- col.idx + 1 # indices (in R) for which variables have scale > 1e-6
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  a <- a[ind]
  b <- b[, ind, drop = FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  
  if (warn & any(iter == max.iter))
    warning("Algorithm failed to converge for some values of lambda")
  
  ## Unstandardize coefficients:
  beta <- Matrix(0, nrow = p, ncol = length(lambda), sparse = TRUE)
  bb <- b / scale[col.idx]
  beta[col.idx, ] <- bb
  aa <- a - as.numeric(crossprod(center[col.idx], bb))
  
  ## Names
  names(aa) <- colnames(beta) <- round(lambda, digits = 4)
  
  ## Output
  structure(list(
    intercept = aa,
    beta = beta,
    iter = iter,
    lambda = lambda,
    penalty = penalty,
    family = "binomial",
    alpha = alpha,
    loss = loss,
    penalty.factor = penalty.factor,
    n = n,
    p = p,
    center = center,
    scale = scale,
    y = yy,
    col.idx = col.idx,
    rejections = rejections
  ), class = "big_sp")
}