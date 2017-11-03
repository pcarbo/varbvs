\name{predict.varbvs}
\alias{predict.varbvs}
\title{Make predictions from a model fitted by varbvs.}
\description{This function predicts outcomes (Y) given the observed
  variables (X) and observed covariates (Z), and a model fitted using
  \code{\link{varbvs}}.}
\usage{
  \method{predict}{varbvs}(object, X, Z = NULL,
                           type = c("link","response","class"),
                           averaged = TRUE, \dots)
}
\arguments{
  \item{object}{Output of function \code{\link{varbvs}}.}

  \item{X}{n x p input matrix, in which p is the number of variables,
    and n is the number of samples for which predictions will be made
    using the fitted model. X cannot be sparse, and cannot have any
    missing values (NA).}

  \item{Z}{n x m covariate data matrix, where m is the number of
    covariates. Do not supply an intercept as a covariate (i.e., a
    column of ones), because an intercept is automatically included in
    the regression model. For no covariates, set \code{Z = NULL}.}

  \item{type}{Type of prediction to output. The default, "link", gives
    the linear predictors for \code{family = "binomial"}, and gives the
    fitted values for \code{family = "gaussian"}. For logistic
    regression (\code{family = "binomial"}), there are two alternative
    predictions: "response" givees the fitted probabilities, and "class"
    produces the maximum-probability outcome (0 or 1).}

  \item{averaged}{When \code{averaged = TRUE}, the predictions are
    computed by averaging over the hyperparameter settings, treating
    \code{object$logw} as (unnormalized) log-marginal
    probabilities. (See \code{\link{varbvs}} for more details about
    averaging.) When \code{averaged = FALSE}, the predictions are
    returned as a matrix when one row for each data sample, and one
    column for each hyperparameter setting.}
  
  \item{\dots}{Other arguments to generic predict function. These
    extra arguments are not used here.}
}
\details{
  Note that the classification probabilities \eqn{Pr(Y = 1 | X, Z,
  \theta)} are not guaranteed to be calibrated under the variational
  approximation.
}  
\value{
  When \code{averaged = TRUE}, the output is a vector containing the
  predicted outcomes for all samples. For \code{family = "binomial"},
  all vector entries are 0 or 1.

  When \code{averaged = FALSE}, the return value is a matrix with one
  row for each sample, and one column for each hyperparameter setting.  
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis}
\bold{7}, 73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{\link{varbvs}}, \code{\link{summary.varbvs}}}
\examples{
  # See help(varbvs) for examples.
}
