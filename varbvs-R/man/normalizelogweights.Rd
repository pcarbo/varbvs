\name{normalizelogweights}
\alias{normalizelogweights}
\title{Compute normalized probabilities.}
\description{Compute normalized probabilities from unnormalized
  log-probabilities.}
\usage{
  normalizelogweights(logw)
}
\arguments{
  \item{logw}{Vector of unnormalized log-probabilities.}
}
\details{
  Guards against underflow or overflow by adjusting the
  log-probabilities so that the largest probability is 1.
}
\value{Normalized probabilities such that the sum is equal to 1.}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis} \bold{7}, 
73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\examples{
  logw <- rnorm(6)
  w    <- normalizelogweights(logw)
}
