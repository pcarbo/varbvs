\name{varbvs.properties}
\alias{varbvs.properties}
\alias{varbvsmix.properties}
\alias{nobs.varbvs}
\alias{case.names.varbvs}
\alias{variable.names.varbvs}
\alias{labels.varbvs}
\alias{coef.varbvs}
\alias{coef.varbvsmix}
\alias{confint.varbvs}
\alias{fitted.varbvs}
\alias{resid.varbvs}
\alias{residuals.varbvs}
\alias{deviance.varbvs}
\title{Accessing Properties of Fitted varbvs Models}
\usage{
\method{nobs}{varbvs}(object, \dots)
\method{case.names}{varbvs}(object, \dots)
\method{variable.names}{varbvs}(object, full = FALSE,
                                include.threshold = 0.01, \dots)
\method{labels}{varbvs}(object, \dots)
\method{coef}{varbvs}(object, \dots)
\method{coef}{varbvsmix}(object, \dots)
\method{confint}{varbvs}(object, parm, level = 0.95, ...)
\method{fitted}{varbvs}(object, \dots)
\method{resid}{varbvs}(object, type = c("deviance","response"), \dots)
\method{residuals}{varbvs}(object, type = c("deviance","response"), \dots)
\method{deviance}{varbvs}(object, \ldots)
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

  \item{parm}{Confidence intervals are computed for these selected
    variables. These may either be specified as numbers (column indices
    of \code{varbvs} input matrix \code{X}) or names (column names of
    \code{X}). If not specified, confidence intervals will be computed
    for the top 5 variables by posterior inclusion probability.
    Confidence intervals are not provided for covariates (columns of
    \code{Z}); see below for details.}

  \item{level}{Size of confidence level.}
  
  \item{type}{Type of residuals to be returned. This argument is only
    relevant for logistic regression models (\code{family =
    "binomial"}). See \code{\link{varbvs}} for more details about the
    two available types of residuals for logistic regression.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
}
\description{
  All these functions are \code{\link{methods}} for class
  \code{"varbvs"} objects.
}
\details{
  The generic accessor functions \code{nobs}, \code{case.names},
  \code{variable.names} and \code{labels} can be used to extract various
  useful properties of the fitted \code{varbvs} model. Method
  \code{labels}, in particular, returns the names of the candidate
  variables (columns of \code{X}) which may be used, for example, to
  plot posterior inclusion probabilities or effect estimates.

  \code{coef} returns a matrix containing the posterior estimates of the
  regression coefficients at each hyperparameter setting, as well as an
  additional column containing "averaged" coefficient estimates.

  \code{confint} returns confidence intervals (also, equivalently in
  this case, "credible intervals") for all selected variables
  \code{parm}. These are \emph{conditional} confidence intervals; that
  is, conditioned on each variable being included in the regression
  model.

  The \code{confint} return value is different from the usual confidence
  interval (e.g., for an \code{\link[stats]{lm}} result) because a
  confidence interval is provided for each hyperparameter setting, as
  well as an additional "averaged" confidence interval. The confidence
  intervals are returned a list, with one list element per selected
  variable, and each list element is a matrix with columns giving lower
  and upper confidence limits for each hyperparameter setting, as well
  as the averaged limits.

  Note that confidence intervals cannot currently be requested for
  covariates (columns of \code{Z}).
  
  \code{fitted} returns a matrix containing the predicted (or "fitted")
  values of the outcome at each hyperparameter setting. For the logistic
  regression model (\code{family = "binomial"}), each matrix entry gives
  the probability that the binary outcome is equal to 1.

  Likewise, \code{resid} and \code{residuals} each return a matrix
  containing the model residuals at each hyperparameter setting.

  \code{deviance} returns the deviance for the fitted model at each
  hyperparameter setting.
}
\seealso{
  \code{\link{varbvs}}, \code{\link[stats]{nobs}},
  \code{\link[stats]{case.names}}, 
  \code{\link[stats]{variable.names}}, \code{\link{labels}},
  \code{\link[stats]{coef}}, \code{\link[stats]{coef}},
  \code{\link[stats]{fitted}}, \code{\link[stats]{residuals}},
  \code{\link[stats]{deviance}}
}
