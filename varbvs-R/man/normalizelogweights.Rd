\name{normalizelogweights}
\alias{normalizelogweights}
\title{Compute normalized importance weights.}
\description{Compute normalized importance weights from unnormalized
  log-importance weights.}
\usage{
  normalizelogweights(logw)
}
\arguments{
  \item{logw}{Vector of unnormalized log-importance weights.}
}
\details{
  Guards against underflow or overflow by adjusting the log-importance
  weights so that the largest importance weight is one.
}
\value{Normalized importance weights such that the sum of the normalized
importance weights is equal to 1.}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. Bayesian Analysis 7(1), pages
73-108.}
\author{Peter Carbonetto <peter.carbonetto@gmail.com>}
\examples{
  logw <- rnorm(6)
  w    <- normalizelogweights(logw)
}
