% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{logit.CMSA}
\alias{logit.CMSA}
\title{Sparse Logistic Regression}
\usage{
logit.CMSA(G, pheno.all, covar.all, ind.train, ind.test, method,
  alphas = c(1, 0.5, 0.05, 0.001), ncores = nb_cores(),
  family = "logistic")
}
\arguments{
\item{G}{A \link[=FBM.code256-class]{FBM.code256}
(typically \code{<bigSNP>$genotypes}).\cr
\strong{You shouldn't have missing values in your data or SNPs with 0 MAF.}}

\item{pheno.all}{All phenotypes.}

\item{covar.all}{All covariables.}

\item{ind.train}{Indices corresponding to the training set.}

\item{ind.test}{Indices corresponding to the test set.}

\item{method}{Name of the method.}

\item{alphas}{See \link[bigstatsr:big_spLogReg]{bigstatsr::big_spLogReg}.}

\item{ncores}{Number of cores to use. Default is \code{bigstatsr::nb_cores()}.}

\item{family}{Default is \code{"logistic"}. Can also use \code{"linear"}.}
}
\value{
A data frame with 5 variables:
\itemize{
\item "method": the name of the method,
\item "pred": the predictions for the test set,
\item "timing": the execution time,
\item "alpha": the alpha that maximized prediction for the validation sets,
\item "set": the set of non-zero coefficients.
}
}
\description{
Sparse Logistic Regression
}
