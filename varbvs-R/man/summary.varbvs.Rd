\name{summary.varbvs,print.summary.varbvs}
\alias{summary.varbvs}
\alias{print.summary.varbvs}
\title{Summarize a fitted variable selection model.}
\description{Generate a summary of the Bayesian variable selection model
  fitted using variational approximation methods.
}
\usage{
  \method{summary}{varbvs}(object, cred.int = 0.95, nv = 5, nr = 1000, \ldots)
  \method{print}{summary.varbvs}(x, digits = 3, \ldots)
}
\arguments{
  \item{object}{Output of function \code{\link{varbvs}}.}
  \item{cred.int}{Size of credible interval (number between 0 and 1).}
  \item{nv}{Show detailed statistics for top nv variables, ranked
  according to their posterior inclusion probabilities.}
  \item{nr}{Number of Monte Carlo samples to draw to estimate credible
    intervals for coefficients of selected variables.}
  \item{x}{Output of function \code{summary.varbvs}.}
  \item{digits}{Number of digits shown when printing posterior
    probabilities of top nv variables.}
  \item{...}{Additional print arguments.}
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
  value averaged over hyperparameter settings, treating \code{objectlogw}
  as (unnormalized) log-marginal probabilities. The next column, labeled
  "Pr>x", where \code{x = cred.int} gives the credible interval based on
  these weights (computed using function \code{\link{cred}}).

  The third part summarizes the variable selection results. This
  includes the total number of variables included in the model at
  different posterior probability thresholds, and a more detailed
  summary of the variables included in the model with highest posterior
  probability. For \code{family = "gaussian"}, the "PVE" column gives
  the estimated proportion of variance in the outcome explained by the
  variable (conditioned on being included in the model). Again, this is
  the PVE after removing variance in the outcome due to linear effects
  of the covariates. Finally, note that the credible intervals reported
  in the right-most column are Monte Carlo estimates, so the interval
  will be slightly different each time \code{summary} is called; a
  more accurate estimate can be obtained by setting input \code{nr} to a
  larger number, at the cost of increased computation time.  
}  
\value{
  An object of class \code{summary.varbvs}, to be printed by
  \code{print.summary.varbvs}.
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis}
\bold{7}, 73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{\link{varbvs}}, \code{\link{varbvscoefcred}}}
\examples{
  # See help(varbvs) for examples.
}
