\name{varbvs.properties}
\alias{nobs.varbvs}
\alias{variable.names.varbvs}
\alias{labels.varbvs}
\alias{coef.varbvs}
\title{Accessing Properties of Fitted varbvs Models}
\usage{
\method{nobs}{varbvs}(object, \dots)
\method{variable.names}{varbvs}(object, full = FALSE,
                                include.threshold = 0.01, \dots)
\method{labels}{varbvs}(object, \dots)
\method{coef}{varbvs}(object, \dots)
}
\arguments{
  \item{object}{An object inheriting from class \code{varbvs}, usually
    the result of calling function \code{\link{varbvs}}.}

  \item{full}{logical; if \code{TRUE}, names of all variables (columns
    of \code{X}) are returned, including variables that have zero
    probability of being included in the regression model.}
  
  \item{include.threshold}{When \code{full = FALSE}, names of all
    variables (columns of \code{X}) with "averaged" posterior inclusion
    probability greater than \code{include.threshold} are returned.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
}
\description{
  All these functions are \code{\link{methods}} for class
  \code{"varbvs"} objects.
}
\details{
  The generic accessor functions \code{nobs}, \code{variable.names},
  \code{labels}, ... can be used to extract various useful properties of
  the fitted \code{varbvs} model.

  Method \code{labels} returns the names of the candidate variables
  (columns of \code{X}) which may be used, for example, to plot
  posterior inclusion probabilities or effect estimates.

  \code{coef} returns a matrix containing the the posterior estimates of
  the regression coefficients at each hyperparameter setting, as well as
  an additional column containing "averaged" coefficient estimates.  
}
\seealso{
  \code{\link{varbvs}}, \code{\link[stats]{nobs}},
  \code{\link[stats]{variable.names}}, \code{\link{labels}},
  \code{\link[stats]{coef}}
}
