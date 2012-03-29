\name{varbvsoptimize}
\alias{varbvsoptimize}
\alias{varbvsupdate}

\title{Coordinate ascent for variational approximation to Bayesian
       variable selection in linear regression}

\description{\code{varbvsoptimize} implements the fully-factorized
variational approximation for Bayesian variable selection in linear
regression. It finds the "best" fully-factorized variational
approximation to the posterior distribution of the coefficients in a
linear regression model of a continuous outcome (quantitiative trait),
with spike and slab priors on the coefficients. By "best", we mean the
approximating distribution that locally minimizes the Kullback-Leibler
divergence between the approximating distribution and the exact
posterior.

\code{varbvsupdate} runs a single iteration of the coordinate ascent
updates to maximize the variational lower bound or, equivalently,
to minimize the Kullback-Leibler divergence objective.}

\usage{

varbvsoptimize(X, y, sigma, sa, logodds, alpha0 = NULL,
mu0 = NULL, verbose = TRUE)

varbvsupdate(X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, S)
}

\arguments{

  \item{X}{Matrix of observations about the variables (or features). It
           has \code{n} rows and \code{p} columns, where \code{n} is the
           number of samples, and \code{p} is the number of variables.}

  \item{sigma}{Scalar giving the variance of the residual.}

  \item{sa}{\code{sa*sigma} is the prior variance of the regression
            coefficients.}

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
  
  \item{y}{Vector of observations about the outcome. It is a vector of
           length \code{n}.}

  \item{verbose}{Set \code{verbose = FALSE} to turn off reporting the 
                 algorithm's progress.}

  \item{xy}{Equal to \code{t(X) \%*\% y}, where \code{y} is the vector of
            observations about the outcome.}

  \item{d}{Equal to \code{diag(t(X) \%*\% X)}.}
	  
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

To account for an intercept, \code{y} and \code{X} must be centered
beforehand so that \code{y} and each column of \code{X} has a mean of
zero.

The computational complexity of running \code{varbvsupdate} is
\code{O(n*length(S))}. For efficient computation, most of the work is
done by \code{varbvsupdateR}, a function implemented in C. The call to
the C function in the shared library is made using \code{\link{.C}}.

}

\value{

\code{varbvsoptimize} returns a list containing four components:
variational parameters \code{alpha}, \code{mu} and \code{s}, and the
variational estimate of the marginal log-likelihood \code{lnZ}.

\code{varbvsupdate} returns a list containing three components:
the updated variational parameters \code{alpha} and \code{mu}, and the
updated matrix-vector product \code{Xr}.}

\author{Peter Carbonetto}

\seealso{
\code{\link{varsimbvs}}
}

\examples{

## Randomly generate genotypes and quantitative trait measurements for
## 1000 SNPs and 500 individuals, in which the variance of the residual
## error is 4. Of these SNPs, 10 have nonzero additive effects on the
## trait.
se   <- 4
snps <- create.snps(1000,10)
data <- create.data(snps$maf,snps$beta,se,500)

## Compute the variational approximation given (appropriate) choices for
## the hyperparameters.
result <- varbvsoptimize(data$X,data$y,4,1/4,log(0.01/0.99))

## View the posterior inclusion probabilities for the (true) causal SNPs.
S <- which(snps$beta != 0)
cbind(snps$beta[S],result$alpha[S])

## View the largest posterior inclusion probability for a SNP that has
## no effect on the quantitative trait.
S <- which(snps$beta == 0)
i <- S[which.max(result$alpha[S])]
cbind(snps$beta[i],result$alpha[i])

}
