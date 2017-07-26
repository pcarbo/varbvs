\name{varbvsindep}
\alias{varbvsindep}
\title{Compute posterior statistics, ignoring correlations.}
\description{Compute the mean and variance of the coefficients, and the
posterior inclusion probabilities (PIPs), ignoring correlations between
variables. This is useful for inspecting or visualizing groups of
correlated variables (e.g., genetic markers in linkage disequilibrium).
}
\usage{
  varbvsindep (fit, X, Z, y)
}
\arguments{
  \item{fit}{Output of function \code{varbvs}.}

  \item{X}{n x p input matrix, where n is the number of samples,
           and p is the number of variables. X cannot be sparse,
           and cannot have any missing values (NA).}

  \item{Z}{n x m covariate data matrix, where m is the number of
           covariates. Do not supply an intercept as a covariate
           (i.e., a column of ones), because an intercept is
           automatically included in the regression model. For no
	   covariates, set \code{Z = NULL}.}

  \item{y}{Vector of length n containing observations of binary
          (\code{family = "binomial"}) or continuous (\code{family =
          "gaussian"}) outcome. For a binary outcome, all entries
          of y must be 0 or 1.}
}
\details{
  For the ith hyperparameter setting, \code{alpha[,i]} is the
  variational estimate of the posterior inclusion probability (PIP) for
  each variable; \code{mu[,i]} is the variational estimate of the
  posterior mean coefficient given that it is included in the model; and
  \code{s[,i]} is the estimated posterior variance of the coefficient
  given that it is included in the model.
}  
\value{
  \item{alpha}{Variational estimates of posterior inclusion
    probabilities for each hyperparameter setting.}

  \item{mu}{Variational estimates of posterior mean coefficients for
    each hyperparameter setting.}

  \item{s}{Variational estimates of posterior variances for each
    hyperparameter setting.}
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis} \bold{7}, 
73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{\link{varbvs}}}
