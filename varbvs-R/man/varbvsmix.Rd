\name{varbvsmix}
\alias{varbvsmix}
\title{Fit linear regression with mixture-of-normals priors using
  variational approximation methods.}
\description{In the original formulation \cite{carbonetto-2012}, each
regression coefficient $\beta_i$ was drawn identically from a
spike-and-slab prior, $\beta_i \sim \pi N(0,\sigma^2 \sigma_a^2) +
(1-\pi) \delta_0$, where $\pi$ is the prior inclusion probability,
$\sigma_a^2$ is the prior variance of the regression coefficients
(scaled by $\sigma^2$), and $\delta_0$ is the delta-mass, or
``spike,'' at zero. Here, we instead formulate the ``slab'' as a
mixture of normals.}
\usage{
  varbvsmix(X, Z, y, sa, sigma, w, alpha, mu, update.sigma, update.sa,
            update.w, w.penalty, tol = 1e-4, maxiter = 1e4, verbose = TRUE)
}
\arguments{
  \item{X}{n x p input matrix, where n is the number of samples, and p
           is the number of variables. X cannot be sparse, and cannot
           have any missing values (NA).}

  \item{Z}{n x m covariate data matrix, where m is the number of
           covariates. Do not supply an intercept as a covariate (i.e.,
           a column of ones), because an intercept is automatically
           included in the regression model. For no covariates, set
           \code{Z = NULL}. The covariates are assigned an improper,
           uniform prior. Although improper priors are generally not
           advisable because they can result in improper posteriors and
           Bayes factors, this choice allows us to easily integrate out
           these covariates.

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
\seealso{\code{\link{varbvs}}}
\examples{
  # Add at least one example here.
}
