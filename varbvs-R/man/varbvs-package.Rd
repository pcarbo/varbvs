\name{varbvs-package}
\alias{varbvs-package}
\alias{varbvs}
\docType{package}
\title{Variational inference for Bayesian variable selection}

\description{The \pkg{varbvs} package implements the variational
inference procedure for Bayesian variable selection, as described in the
"Scalable variational inference for Bayesian variable selection in
regression, and its accuracy in genetic association studies" (Bayesian
Analysis 7, March 2012, pages 73-108). This software has been used to
implement Bayesian variable selection for large problems with over a
million variables and thousands of samples.}

\details{

The main functions that implement the variational inference for Bayesian
variable selection are \code{\link{varbvsoptimize}} (for linear
regression) and \code{\link{varbvsbinoptimize}} (for logistic
regression). Function \code{\link{varsimbvs}} demonstrates how the
variational approximation can be combined with importance sampling to
compute the posterior distributions for coefficients and
hyperparameters. Accompanying functions that might be useful to users
include \code{\link{grid3d}}, \code{\link{create.snps}} and
\code{\link{create.data}}.
  
To see a full example of how the variational inference procedure can be
used to compute posterior inclusion probabilities and posterior
distributions of the hyperparameters in a linear regression setting, run
\code{demo(example1)}.
  
To see more details about the \pkg{varbvs} package, including the
license and a list of available functions, use
\code{help(package="varbvs")}.

}

\author{Peter Carbonetto <pcarbo@uchicago.edu>}

\references{Peter Carbonetto and Matthew Stephens (2012). Scalable
variational inference for Bayesian variable selection in regression, and
its accuracy in genetic association studies. Bayesian Analysis 7(1),
pages 73-108.}


