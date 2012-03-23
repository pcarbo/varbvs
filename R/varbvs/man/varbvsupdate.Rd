  #
  #

\name{varbvsupdate}
\alias{varbvsupdate}
\title{Coordinate ascent updates for variational approximation to
       Bayesian variable selection in linear regression}

\description{Runs a single iteration of the coordinate ascent updates to
maximize the variational lower bound for Bayesian variable selection in
linear regression. It adjusts the fully-factorized variational
approximation to the posterior distribution of the coefficients in a
linear regression model of a continuous outcome (quantitiative trait),
with spike and slab priors on the coefficients.}

\usage{varbvsupdate(X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, S)}

\arguments{

  \item{X}{Matrix of observations about the variables (or features). It
           has \code{n} rows and \code{p} columns, where \code{n} is the
           number of samples, and \code{p} is the number of
           variables. Must be in double-precision.}

  \item{sigma}{Scalar giving the variance of the residual.}

  \item{sa}{\code{sa*sigma} is the prior variance of the regression
            coefficients.}

  \item{logodds}{Prior log-odds of inclusion for each variable. It is
                 equal to \code{logodds = log(q/(1-q))}, where \code{q}
                 is the prior probability that each variable is included
                 in the linear model of Y. It is a vector of length
                 \code{p}.}
	    
  \item{xy}{Equal to \code{t(X) \%*\% y}, where \code{y} is the vector of
            observations about the outcome.}

  \item{d}{Equal to \code{diag(t(X) \%*\% X)}.}
	  
  \item{alpha0}{Current variational estimates of the posterior inclusion
                probabilities. It is a vector of length \code{p}.}
  
  \item{mu0}{Current variational estimates of the posterior mean
             coefficients. It is a vector of length \code{p}.}
  
  \item{Xr0}{Equal to \code{X \%*\% (alpha0*mu0)}.}
  
  \item{S}{Order in which the coordinates are updated. It is a vector of
           any length. Each entry of \code{S} must be an integer between
           1 and \code{p}.}}
       
\details{

To account for an intercept, \code{y} and \code{X} must be centered
beforehand so that \code{y} and each column of \code{X} has a mean of
zero.

The computational complexity of this procedure is \code{O(n*length(S))}.
}

\value{
Returns a list containing three components:
\item{alpha}{updated variational estimates of the posterior inclusion
             probabilities.}
\item{mu}{updated variational estimates of the posterior mean
          coefficients.}
\item{Xr}{equal to \code{X \%*\% (alpha*mu)}.}
}

\author{Peter Carbonetto}

\seealso{\code{\link{varbvsoptimize}}}

