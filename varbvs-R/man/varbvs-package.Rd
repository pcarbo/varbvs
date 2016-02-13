\name{varbvs-package}
\alias{varbvs-package}
\docType{package}
\title{Large-scale Bayesian variable selection using variational methods}
\description{Fast algorithms for fitting Bayesian variable selection
  models and computing Bayes factors, in which the outcome (or response
  variable) is modeled using a linear regression or a logistic
  regression. The algorithms are based on the variational approximations
  described in "Scalable variational inference for Bayesian variable
  selection in regression, and its accuracy in genetic association
  studies" (P. Carbonetto and M. Stephens, Bayesian Analysis 7, 2012,
  pages 73-108). This software has been applied to large data sets with
  over a million variables and thousands of samples.}
\details{

The main functionality of this package is implemented in function
\code{\link{varbvs}}. This function selects the most appropriate
algorithm for the data set and selected model (linear or logistic
regression). See \code{help(varbvs)} for details. The varbvs interface
is intended to resemble interface for \pkg{glmnet}, the popular package
for fitting genealized linear models.
  
For more details about the this package, including the license and a
list of available functions, see \code{help(package=varbvs)}.

}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis} \bold{7}, 
73--108.}
