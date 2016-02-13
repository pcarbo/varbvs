\name{varbvscoefcred}
\alias{varbvscoefcred}
\title{Compute credible intervals for regression coefficients.}
\description{Compute Monte Carlo estimates of credible intervals for
coefficients in the fitted variable selection model. This function is
used by summary.varbvs to generate credible intervals for coefficients
of top-ranked variables.
}
\usage{
  varbvscoefcred (fit, vars, cred.int = 0.95, nr = 1000)
}
\arguments{
  \item{fit}{Output of function \code{\link{varbvs}}.}

  \item{vars}{Vector of indices or names of variables. If not specified,
    credible intervals are computed for all variables.}

  \item{cred.int}{Size of credible interval (number between 0 and 1).}

  \item{nr}{Number of Monte Carlo samples to draw to estimate credible
    intervals. Amore accurate estimate of the credible interval can be
    obtained by setting nr to a larger number, at the cost of increased
    computation time.}
}
\details{
  Here, the credible interval [a,b] is simply defined as \code{a =
  quantile(x,0.5 - cred.int/2)} and \code{b = quantile(x,0.5 +
  cred.int/2)}, in which \code{x} is a vector of samples drawn from the
  posterior distribution. 
}  
\value{
  \item{a}{Credible interval lower bounds.}
  \item{b}{Credible interval upper bounds.}
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis} \bold{7},
73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{varbvs}, \code{summary.varbvs}}
