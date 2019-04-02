\name{predict.varbvsmix}
\alias{predict.varbvsmix}
\title{Make predictions from a model fitted by varbvsmix.}
\description{This function predicts outcomes (Y) given the observed
  variables (X) and observed covariates (Z), and a model fitted using
  \code{\link{varbvsmix}}.}
\usage{
  \method{predict}{varbvsmix}(object, X, Z = NULL, \dots)
}
\arguments{
  \item{object}{Output of function \code{\link{varbvsmix}}.}

  \item{X}{n x p input matrix, in which p is the number of variables,
    and n is the number of samples for which predictions will be made
    using the fitted model. X cannot be sparse, and cannot have any
    missing values (NA).}

  \item{Z}{n x m covariate data matrix, where m is the number of
    covariates. Do not supply an intercept as a covariate (i.e., a
    column of ones), because an intercept is automatically included in
    the regression model. For no covariates, set \code{Z = NULL}.}
  
  \item{\dots}{Other arguments to generic predict function. These
    extra arguments are not used here.}
}
\seealso{\code{\link{varbvsmix}}}
\examples{
  # See help(varbvsmix) for examples.
}
