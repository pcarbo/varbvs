\name{varbvsbinoptimize}
\alias{varbvsbinoptimize}
\alias{varbvsbinupdate}

\title{Coordinate ascent for variational approximation to Bayesian
       variable selection in logistic regression}

\description{\code{varbvsbinoptimize} implements the fully-factorized
variational approximation for Bayesian variable selection in logistic
regression. It finds the "best" fully-factorized variational
approximation to the posterior distribution of the coefficients in a
logistic regression model of a binary outcome (e.g. disease status in a
case-control study), with spike and slab priors on the coefficients. By
"best", we mean the approximating distribution that locally minimizes
the Kullback-Leibler divergence between the approximating distribution
and the exact posterior.

\code{varbvsbinupdate} runs a single iteration of the coordinate ascent
updates to maximize the variational lower bound or, equivalently,
to minimize the Kullback-Leibler divergence objective.}

\usage{

varbvsbinoptimize(X, y, sa, logodds, alpha0 = NULL, mu0 = NULL,
                  eta0 = NULL, fixed.eta = FALSE, verbose = TRUE)

varbvsbinupdate(X, sa, logodds, stats, alpha0, mu0, Xr0, S)
}

\arguments{
  
  \item{X}{Matrix of observations about the variables (or features). It
           has \code{n} rows and \code{p} columns, where \code{n} is the
           number of samples, and \code{p} is the number of variables.}

  \item{sa}{Prior variance of the regression coefficients.}

  \item{logodds}{Prior log-odds of inclusion for each variable. It is
                 equal to \code{logodds = log(q/(1-q))}, where \code{q}
                 is the prior probability that each variable is included
                 in the linear model of Y. It may either be a scalar, in
                 which case all the variables have the same prior
                 inclusion probability, or it may be a vector of length
                 \code{p}.}

  \item{alpha0}{Initial variational estimate of posterior inclusion 
                probabilities. It is a vector of length \code{p}. If
                \code{alpha0 = NULL}, the variational parameters are
                initialized at random.}
	      
  \item{mu0}{Initial variational estimate of posterior mean
                coefficients. It is a vector of length \code{p}. If
                \code{mu0 = NULL}, the variational parameters are
                randomly initialized.}
	      
  \item{y}{Vector of observations about the binary outcome. It is a
           vector of length \code{n}.}
  
  \item{eta0}{Initial free parameters specifying the variational lower
              bound on the logistic regression factors. If \code{eta0 =
              NULL}, the free parameters are initialized to a vector of
              ones.}
	      
  \item{fixed.eta}{Set \code{fixed.eta = TRUE} to prevent ETA from being
                   updated.}

  \item{verbose}{Set \code{verbose = FALSE} to turn off reporting the 
                 algorithm's progress.}

  \item{stats}{This is the return value of
  \code{updatestats(X,y,eta)}. See \code{\link{updatestats}} for more
  information on this argument.}

  \item{Xr0}{Equal to \code{X \%*\% (alpha0*mu0)}.}
  
  \item{S}{Order in which the coordinates are updated. It is a vector of
           any length. Each entry of \code{S} must be an integer between
           1 and \code{p}.}
}
       
\details{

\code{alpha}, \code{mu} and \code{s} are the parameters of the
variational approximation and, equivalently, variational estimates of
posterior quantites: under the variational approximation, the ith
regression coefficient is normal with probability \code{alpha[i]};
\code{mu[i]} and \code{s[i]} are the mean and variance of the
coefficient given that it is included in the model. \code{alpha},
\code{mu} and \code{s} are always column vectors of length \code{p}.

Unlike \code{\link{varbvsoptimize}}, \code{y} and \code{X} must not be
centered. Instead, we will account for the intercept as we update the
variational approximation.

Also note that the residual variance parameter \code{sigma} is not
needed to model a binary trait.

The computational complexity of running \code{varbvsbinupdate} is
\code{O(n*length(S))}. For efficient computation, most of the work is
done by \code{varbvsbinupdateR}, a function implemented in C. The call to
the C function in the shared library is made using \code{\link{.C}}.  }

\value{

\code{varbvsbinoptimize} returns a list containing five components:
variational parameters \code{alpha}, \code{mu} and \code{s}, the vector
of free parameters \code{eta}, and the variational estimate of the
marginal log-likelihood \code{lnZ}.

\code{varbvsbinupdate} returns a list containing three components:
the updated variational parameters \code{alpha} and \code{mu}, and the
updated matrix-vector product \code{Xr}.}

\author{Peter Carbonetto}

