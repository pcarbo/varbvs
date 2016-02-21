\name{pve2sa}
\alias{pve2sa}
\title{Transform PVE estimates to prior variance estimates.}
\description{
PVE2SA(SX,H,THETA0) returns the prior variance of the additive effects on
the outcome given the prior estimate of the proportion of variance
explained by the additive effects (H). THETA0 is the (base 10) logarithm
of the prior odds for inclusion, and SX is the sum of the sample variances
for all the explanatory variables.
NOTE: Only for linear regression model.
}
\usage{
  pve2sa (X, h, logodds)
}
\arguments{
  \item{X}{n x p input matrix, where n is the number of samples,
           and p is the number of variables. X cannot be sparse,
           and cannot have any missing values (NA).}

  \item{h}{Proportion of variance in continuously valued outcome that is
  explained by variables. May be a single number or a vector.}

  \item{logodds}{Prior log-odds that a variable is included in the
    regression model. May be a single number of a vector with the same
    number of entries as \code{h}. See \code{\link{varbvs}} for more
    details on this parameter}.
}
\value{
  The return value is a vector with as many entries as inputs \code{h}
  and \code{logodds}.
}
\references{
  Y. Guan and M. Stephens (2011). Bayesian variable selection regression
  for genome-wide association studies and other large-scale
  problems. \emph{Annals of Applied Statistics} \bold{5}, 1780–-1815.

  X. Zhou, P. Carbonetto and M. Stephens (2013). Polygenic modeling with
  Bayesian sparse linear mixed models. \emph{PLoS Genetics} \bold{9},
  e1003264.
}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{varbvs}}
