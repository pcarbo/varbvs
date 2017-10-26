\name{varbvsproxybf}
\alias{varbvsproxybf}
\title{Compute Bayes factors measuring improvement-in-fit along 1
  dimension.}
\description{For each candidate variable j, this function returns a
  Bayes factor measuring the improvement in fit when variable j is
  included in the model instead of variable i; that is, a larger Bayes
  factor indicates a better model fit by swapping variables i and
  j. From an optimization perspective, this could be viewed as
  addressing the following question: if you had to update the
  variational parameters for one variable so as to improve the "fit" of
  the variational approximation after setting the posterior inclusion
  probability for variable i to zero, which variable would you choose?
}
\usage{
  varbvsproxybf(X, Z, y, fit, i, vars)
}
\arguments{
  \item{X}{n x p input matrix, where n is the number of samples, and p
    is the number of variables. X cannot be sparse, and cannot have any
    missing values (NA).}

  \item{Z}{n x m covariate data matrix, where m is the number of
    covariates. Do not supply an intercept as a covariate (i.e., a
    column of ones), because an intercept is automatically included in
    the regression model. For no covariates, set \code{Z = NULL}.}

  \item{y}{Vector of length n containing values of the continuous
    outcome.}

  \item{fit}{An object inheriting from class \code{varbvs}, usually the
    result of calling function \code{\link{varbvs}}. Currently, this is
    only implemented for linear regression (\code{family =
    "gaussian"}); any other choice will produce an error.}

  \item{i}{Variable against will. Typically, will be a variable included
    in the regression model with high probability, but not always.}

  \item{vars}{Set of candidate "proxy" variables. This set may include
    \code{i}, but not does not have to.}
}
\value{
  
  \code{varbvsproxybf} returns a list with the following components:

  \item{BF}{Matrix containing Bayes factors for each candidate proxy
    variable and for each hyperparameter setting.}

  \item{mu}{Matrix containing estimated posterior means for each
    candidate proxy variable and for each hyperparameter setting.}

  \item{s}{Matrix containing estimated posterior variances for each
    candidate proxy variance for each hyperparameter setting.}
}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{\link{varbvs}}}
