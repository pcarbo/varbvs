\name{rand,randn}
\alias{rand}
\alias{randn}
\title{Return matrices of pseudorandom values.}
\description{Generate matrices of pseudorandom values.}
\usage{
  rand(m,n)
  randn(m,n)
}
\arguments{
  \item{m}{Number of matrix rows.}
  \item{n}{Number of matrix columns.}
}
\details{
  Function \code{rand} returns a matrix containing pseudorandom values
  drawn from the standard uniform distribution (using \code{runif}).
  Function \code{randn} returns a matrix containing pseudorandom values
  drawn from the standard normal distribution (using \code{rnorm}).
}
\value{An m x n numeric matrix.}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis} \bold{7}, 
73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\examples{
  x <- rand(10,5)
  y <- randn(10,5)
}
