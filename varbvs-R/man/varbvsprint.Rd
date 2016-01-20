\name{varbvsprint}
\alias{varbvsprint}
\title{Print a fitted variable selection model.}
\description{Print a summary of the Bayesian variable selection model
  fitted using variational approximation methods.
}
\usage{
  varbvsprint(fit, cred.int = 0.95, n = 5, nr = 1000)
}
\arguments{
  \item{fit}{Output of function \code{varbvs}.}

  \item{cred.int}{Size of credible interval (number between 0 and 1).}

  \item{n}{Show detailed statistics for top n > 0 variables, ranked
  according to their posterior inclusion probabilities.}

  \item{nr}{Number of Monte Carlo samples to draw to estimate the
                proportion of variance in the outcome explained by the
                fitted model.}
}
\details{

  The printed summary is divided into three parts. The first part
  summarizes the data and optimization settings. It also reports the
  hyperparameter setting that yields the largest marginal
  likelihood---more precisely, the approximate marginal likelihood
  computed using the variational method. For the linear regression only
  (\code{family = "gaussian"}), it reports the estimated proportion of
  variance in the outcome explained by the model (PVE), and the credible
  interval of the PVE estimate brackets. Note that this is the PVE in
  the outcome after removing variance in the outcome due to linear
  effects of the covariates.

  The second part summarizes the approximate posterior distribution of
  the hyperparameters (sigma, sa, logodds). The "estimate" column is the
  value averaged over hyperparameter settings, treating \code{fit$logw}
  as (unnormalized) log-importance weights. The next column, labeled
  "Pr>x", where \code{x = cred.int} gives the credible interval based on
  these weights (computed using function \code{cred}).

  The third part summarizes the variable selection results. This
  includes the total number of variables included in the model at
  different posterior probability thresholds, and a more detailed
  summary of the variables included in the model with highest posterior
  probability. For \code{family = "gaussian"}, the "PVE" column gives
  the estimated proportion of variance in the outcome explained by the
  variable (conditioned on being included in the model). Again, this is
  the PVE after removing variance in the outcome due to linear effects
  of the covariates.
  
}  
\value{
  Input \code{fit} is silently returned.
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. Bayesian Analysis 7(1), pages
73-108.}
\author{Peter Carbonetto <peter.carbonetto@gmail.com>}
\seealso{\code{varbvsprint}, \code{varbvscoefcred}, \code{varbvspve},
  \code{varbvsnorm}, \code{varbvsbin}, \code{varbvsbinz},
  \code{varbvspve} and \code{normalizelogweights}
}
\examples{
See help(varbvs) for examples.
}
