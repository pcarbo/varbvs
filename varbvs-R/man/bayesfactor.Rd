\name{bayesfactor}
\alias{bayesfactor}
\title{Compute numerical estimate of Bayes factor.}
\description{
Computes numerical estimate of
\deqn{
  BF = Pr(data | H1) / Pr(data | H0),
}
the probability of the data given the "alternative" hypothesis (H1) over
the probability of the data given the "null" hypothesis (H0). This is
also known as a Bayes factor (see Kass & Raftery, 1995). Here we assume
that although these probabilities cannot be computed analytically
because they involve intractable integrals, we can obtain reasonable
estimates of these probabilities with a simple numerical approximation
over some latent variable assuming the prior over this latent variable
is uniform. The inputs are the log-probabilities
\deqn{
  Pr(data, Z0 | H0) = Pr(data | Z0, H0) x Pr(Z0 | H0),
  Pr(data, Z1 | H1) = Pr(data | Z1, H1) x Pr(Z1 | H1),
}
where Pr(Z0 | H0) and Pr(Z1 | H1) are uniform over all Z0 and Z1.

Alternatively, this function can be viewed as computing an importance
sampling estimate of the Bayes factor; see, for example, R. M. Neal,
"Annealed importance sampling", Statistics and Computing, 2001. This
formulation described above is a special case of importance sampling
when the settings of the latent variable Z0 and A1 are drawn from the
same (uniform) distribution as the prior, Pr(Z0 | H0) and Pr(Z1 | H1),
respectively.
}
\usage{
  bayesfactor (logw0, logw1)
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
  genetic association studies. \emph{Bayesian Analysis} \bold{7},
  73--108.

  R. E. Kass and A. E. Raftery (1995). Bayes Factors. \emph{Journal of the
  American Statistical Association} \bold{90}, 773--795.

  R. M. Neal (2001). Annealed importance sampling. \emph{Statistics and
  Computing} \bold{11}, 125--139.
}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{varbvs}, \code{normalizelogweights}}
