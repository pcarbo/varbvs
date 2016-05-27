\name{varbvs-internal}
\title{Internal varbvs functions}
\alias{logit}
\alias{logpexp}
\alias{sigmoid}
\alias{logsigmoid}
\alias{diagsq}
\alias{diagsqt}
\alias{updatestats_varbvsbin}
\alias{updatestats_varbvsbinz}
\alias{update_eta}
\alias{update_etaz}
\alias{varbvsnorm}
\alias{varbvsbin}
\alias{varbvsbinz}
\alias{varbvsnormupdate}
\alias{varbvsbinupdate}
\alias{varbvsbinzupdate}
\alias{varbvspve}
\alias{var1}
\alias{var1.cols}
\description{Internal varbvs functions}
\usage{
diagsq(X,a = NULL)
diagsqt(X,a = NULL)
var1(x)
var1.cols(X)
varbvspve(X,fit,nr = 1000)
varbvsnorm(X,y,sigma,sa,logodds,alpha,mu,tol = 1e-4,maxiter = 1e4,
           verbose = TRUE,outer.iter = NULL,update.sigma = TRUE,
           update.sa = TRUE,n0 = 0,sa0 = 0)
varbvsbin(X,y,sa,logodds,alpha,mu,eta,tol = 1e-4,maxiter = 1e4,
          verbose = TRUE,outer.iter = NULL,update.sa = TRUE,
          optimize.eta = TRUE,n0 = 0,sa0 = 0)
varbvsbinz(X,Z,y,sa,logodds,alpha,mu,eta,tol = 1e-4,maxiter = 1e4,
           verbose = TRUE,outer.iter = NULL,update.sa = TRUE,
           optimize.eta = TRUE,n0 = 0,sa0 = 0)
}
\author{Peter Carbonetto <peter.carbonetto@gmail.com>}

\details{These functions are only intended to be used by expert
  users. Here we provide brief descriptions of some of these internal
  functions. \code{diagsq(X,a)} and \code{diagsqt(X,a)} compute the
  diagonal entries of the matrix products \eqn{X'*A*X} and \eqn{X*A*X'},
  respectively, where \code{A = diag(a)}, but the computation is done
  more efficiently, and without having to store an intermediate matrix
  of the same size as the matrix.

  \code{var1(x)} returns the second moment of vector x about its mean.
  
  \code{var1.cols(X)} computes the second moment of each column of X about
  its mean.

  \code{varbvspve} draws posterior estimates of the proportion of
  variance in Y explained by the Bayesian variable selection model
  fitted using a variational approximation. This function is only valid
  for the linear regression model \code{(family = "gaussian")}.

  Functions \code{varbvsnorm}, \code{varbvsbin} and \code{varbvsbinz}
  implement the co-ordinate ascent algorithm to fit the fully-factorized
  variational approximation for Bayesian variable selection, conditioned
  on settings of the hyperparameters. These functions implement the
  algorithm for the linear regression, logistic regression with an
  intercept, and logistic regression with arbitrary covariates,
  respectively.
}
\keyword{internal}
