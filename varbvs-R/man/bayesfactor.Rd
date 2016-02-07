\name{bayesfactor}
\alias{bayesfactor}
\title{Compute numerical estimate of Bayes factor.}
\description{
Computes numerical estimate of \eqn{BF = Pr(data | H1) / Pr(data | H0),}
the probability of the data given the alternative hypothesis (H1) over
the probability of the data given the null hypothesis (H0). This is also
known as a Bayes factor (see Kass & Raftery, 1995). Here we assume that
although these probabilities cannot be computed analytically because
they involve intractable integrals, we can obtain reasonable estimates
of these probabilities with a simple numerical approximation over some
latent variable, Z, assuming the prior over Z is uniform. The inputs are
the log-probabilities
\eqn{Pr(data, Z | H) = Pr(data | Z, H) x Pr(Z | H),}
where Pr(Z | H) is uniform over all Z, and H is H0 or H1.

Alternatively, this function can be used to compute an importance
sampling estimate of the Bayes factor; see, for example, Neal, 2001.
The importance sampling estimate is equivalent to the simple numerical
approximation when the settings of the latent variable Z are drawn from
the same distribution as the prior, Pr(Z | H).
}
\usage{
  varbvsbayesfactor (logw0, logw1)
}
\arguments{
  \item{logw0}{log-probabilities or log-importance weights under H0.}
  \item{logw1}{log-probabilities or log-importance weights under H1.}
}
\value{
  The estimated Bayes factor.
}
\references{
  P. Carbonetto and M. Stephens (2012). Scalable variational inference
  for Bayesian variable selection in regression, and its accuracy in
  genetic association studies. Bayesian Analysis 7(1), pages 73-108.\cr

  R. E. Kass and A. E. Raftery (1995). Bayes Factors. Journal of the
  American Statistical Association 90(430), pages 773–795.\cr

  R. M. Neal (2001). Annealed importance sampling. Statistics and
  Computing 11(2), pages 125–139.
}
\author{Peter Carbonetto <peter.carbonetto@gmail.com>}
\seealso{\code{normalizelogweights}}
