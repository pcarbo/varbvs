\name{varbvsmix}
\alias{varbvsmix}
\title{Fit linear regression with mixture-of-normals priors using
  variational approximation methods.}
\description{Description goes here.}
\usage{
  varbvsmix(X, Z, y, sa, sigma, w, alpha, mu, update.sigma, update.sa,
            update.w, w.penalty, tol = 1e-4, maxiter = 1e4, verbose = TRUE)
}
\arguments{
  \item{X}{Description of input argument goes here.}

  \item{Z}{Description of input argument goes here.}

  \item{y}{Description of input argument goes here.}

  \item{sa}{Description of input argument goes here.}

  \item{sigma}{Description of input argument goes here.}

  \item{w}{Description of input argument goes here.}

  \item{alpha}{Description of input argument goes here.}

  \item{mu}{Description of input argument goes here.}

  \item{update.sigma}{Description of input argument goes here.}

  \item{update.sa}{Description of input argument goes here.}

  \item{update.w}{Description of input argument goes here.}

  \item{w.penalty}{Description of input argument goes here.}

  \item{tol}{Description of input argument goes here.}

  \item{maxiter}{Description of input argument goes here.}

  \item{verbose}{Description of input argument goes here.}
}
\details{
  Details go here.
}
\value{
  An object with S3 class \code{c("varbvsmix","list")}.

  \item{n}{Number of data samples used to fit model.}

  \item{ncov}{Number of covariates (columns of Z).}
}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\examples{
  # Add at least one example here.
}
