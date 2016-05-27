\name{predict.varbvs}
\alias{predict.varbvs}
\title{Make predictions from a model fitted by varbvs.}

\description{This function predicts outcomes (Y) given the observed
  variables (X) and observed covariates (Z), and a model fitted using
  \code{\link{varbvs}}.}

\usage{
  \method{predict}{varbvs}(object, X, Z = NULL, \ldots)
}
\arguments{
  \item{object}{Output of function \code{\link{varbvs}}.}

  \item{X}{n x p input matrix, in which p is the number of variables,
           and n is the number of samples for which predictions will be
           made using the fitted model. X cannot be sparse, and cannot
           have any missing values (NA).}

  \item{Z}{n x m covariate data matrix, where m is the number of
           covariates. Do not supply an intercept as a covariate
           (i.e., a column of ones), because an intercept is
           automatically included in the regression model. For no
	   covariates, set \code{Z = NULL}.}

  \item{\dots}{Other arguments to generic predict function. These
    extra arguments are not used here.}
}
\details{
  For the logistic regression model, we do not provide classification
  probabilities \eqn{Pr(Y = 1 | X, Z)} because these probabilities are not
  necessarily calibrated under the variational approximation.
  
  The predictions are computed by averaging over the hyperparameter
  settings, treating \code{object$logw} as (unnormalized) log-marginal
  probabilities. See \code{\link{varbvs}} for more details about
  correctly using \code{object$logw} for approximate numerical integration
  over the hyperparameters, for example by treating these as importance
  weights.
}  
\value{
  Vector containing the predicted outcomes for all samples. For
  \code{family = "binomial"}, all vector entries are 0 or 1.
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis}
\bold{7}, 73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{varbvs}, \code{summary.varbvs}}
\examples{
  # See help(varbvs) for examples.
}
