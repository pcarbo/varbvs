\name{cred}
\alias{cred}
\title{Estimate credible interval.}
\description{
  Estimate credible interval from weighted samples.
}
\usage{
  cred(x, x0, w = NULL, cred.int = 0.95)
}
\arguments{
  \item{x}{Vector of random samples of variable.}
  \item{x0}{Mean of median of variable.}
  \item{w}{Weight > 0 assigned to each sample. If w = NULL, all weights
    are the same.}
  \item{cred.int}{Credible interval must contain probability mass of at
    least this amount. A number between 0 and 1.}
}
\details{
Credible interval [a,b] is defined as smallest interval containing x0
that contains \code{cred.int} of the probability mass. Note that the credible
interval is not necessarily symmetric about x0. (Other definitions of the
credible interval are possible.)

The algorithm is quadratic in the length of x (and w), so should not be
used for large numbers of samples.
}
\value{list(a = a,b = b).}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis} \bold{7}, 
73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\examples{
  x   <- rnorm(100)
  out <- cred(x,mean(x),cred.int = 0.68)
}
