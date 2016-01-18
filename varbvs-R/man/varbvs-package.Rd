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
variable selection are *removed* (for linear
regression) and *removed* (for logistic
regression). Function *removed* demonstrates how the
variational approximation can be combined with importance sampling to
compute the posterior distributions for coefficients and
hyperparameters. 
  
To see a full example of how the variational inference procedure can be
used to compute posterior inclusion probabilities and posterior
distributions of the hyperparameters in a linear regression setting, run
*removed*.
  
To see more details about the \pkg{varbvs} package, including the
license and a list of available functions, use
\code{help(package="varbvs")}.
}

\author{Peter Carbonetto <peter.carbonetto@gmail.com>}

\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. Bayesian Analysis 7(1), pages
73-108.}
