\name{varsimbvs}
\alias{varsimbvs}

\title{Full inference procedure for Bayesian variable selection
in linear regression}

\description{\code{varsimbvs} is a special implementation of the
variational inference procedure used in the two simulation studies for
the Bayesian Analysis paper. The main distinguishing feature of this
procedure is the choice of priors for the hyperparameters of the
variable selection model. In addition, we also avoid erratic behaviour
in the variational approximation by first searching for a good
initialization of the variational parameters. This inference procedure
involves an inner loop and an outer loop. The inner loop consists of
running a coordinate ascent algorithm to tighten the variational lower
bound given a setting of the hyperparameters (this inner loop is
implemented by \code{\link{varbvsoptimize}}). The outer loop computes
importance weights for all combinations of the hyperparameters.}

\usage{varsimbvs(X, y, sigma, sa, log10q, a, b, ca)}

\arguments{

  \item{X}{Matrix of observations about the variables (or features). It
           has \code{n} rows and \code{p} columns, where \code{n} is the
           number of samples, and \code{p} is the number of variables.}
	 
  \item{y}{Vector of observations about the outcome. It is a vector of
           length \code{n}.}

  \item{sigma, sa, log10q}{These arguments specify the combinations of
    the hyperparameter settings. These inputs must be vectors of the
    same length. For each combination of the hyperparameters, we compute
    an importance weight. \code{sigma} is the residual variance,
    \code{sigma*sa} is the prior variance of the regression
    coefficients, and \code{log10q} is the (base 10) logarithm of the
    prior inclusion probability.}

  \item{a, b, ca}{These arguments are all positive scalars specifying
    the priors on the hyperparameters. \code{a} and \code{b} are the
    prior sample sizes for the beta prior on the prior inclusion
    probability. We assume a uniform prior on the "proportion of
    variance explained", except that we replace the prior inclusion
    probability in the proportion of variance explained by a constant,
    \code{ca}. This is done purely for convenience, so that
    hyperparameter \code{sa} does not depend on the prior inclusion
    probability a priori, making it easier to implement the Markov chain
    Monte Carlo (MCMC) method (see the Bayesian Analysis article for
    details). We assume the standard noninformative prior on the
    residual variance \code{sigma}.}

}

\value{
Returns a list containing three components:

\item{alpha}{variational estimates of the posterior inclusion
             probabilities, averaged over settings of the
             hyperparameters.}

\item{mu}{variational estimates of the posterior mean coefficients,
          averaged over settings of the hyperparameters.}

\item{w}{normalized importance weights for all settings of the
         hyperparameters}
}
     
\author{Peter Carbonetto}

\seealso{
  \code{\link{varbvsoptimize}}
}

\examples{

## Run the demonstration R Script.
\dontrun{demo(example1)}

}