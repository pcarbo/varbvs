\name{varbvsmix}
\alias{varbvsmix}
\title{Fit linear regression with mixture-of-normals priors using
  variational approximation methods.}
\description{Find the "best" fully-factorized approximation to the
  posterior distribution of the coefficients, with linear regression
  likelihood and mixture-of-normals priors on the coefficients. By
  "best", we mean the approximating distribution that locally minimizes
  the Kullback-Leibler divergence between the approximating distribution
  and the exact posterior. In the original formulation (see
  \code{\link{varbvs}}), each regression coefficient was drawn
  identically from a spike-and-slab prior. Here, we instead formulate
  the ``slab'' as a mixture of normals.}
\usage{
  varbvsmix(X, Z, y, sa, sigma, w, alpha, mu, update.sigma, update.sa,
            update.w, w.penalty, drop.threshold = 1e-8, tol = 1e-4,
            maxiter = 1e4, update.order = 1:ncol(X), verbose = TRUE)
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

  \item{sa}{Vector specifying the prior variance of the regression
    coefficients (scaled by \code{sigma}) for each mixture component. If
    not specified, a "reasonable" set of prior variances is
    automatically selected based on a simple regression analysis in
    which all variables (columns of \code{X}) are treated as
    independent. To override this default, set \code{sa} to vector of
    variances in which the first mixture component is the "spike", and
    therefore should be zero. Alternatively, \code{sa} may be an integer
    greater than 1, in which case this controls the number of mixture
    components, and the variances of the individual mixture components
    are automatically selected as described.}

  \item{sigma}{Residual variance parameter. If missing, it is
    automatically fitted to the data by computing an approximate
    maximum-likelihood estimate.}

  \item{w}{If missing, it is automatically fitted to the data by
    computing an approximate maximum-likelihood estimate.}

  \item{alpha}{Initial estimates of the approximate posterior mixture
    assignment probabilities. These should be specified as a p x K
    matrix, where K is the number of mixture components. Each row must
    add up to 1.}

  \item{mu}{Initial estimates of the approximate regression coefficients
    conditioned on being drawn from each of the K mixture
    components. These estimates should be provided as a p x K matrix,
    where K is the number of mixture components.}

  \item{update.sigma}{If \code{TRUE}, sigma is fitted to data using an
    approximate EM algorithm, in which case argument \code{sigma}, if
    provided, is the initial estimate.}

  \item{update.sa}{Currently, estimate of mixture component variances is
    not implemented, so this must be set to \code{TRUE}, otherwise an
    error will be generated.}

  \item{update.w}{If \code{TRUE}, mixture weights are fitted using an
    approximate EM algorithm, in which case argument \code{w}, if
    provided, is the initial estimate.}

  \item{w.penalty}{Penalty term for the mixture weights. It is useful
    for "regularizing" the estimate of \code{w} when we do not have a
    lot of information. It should be a vector with one positive entry
    for each mixture component. Larger values place more weight on the
    corresponding mixture components. It is based on the Dirichlet
    distribution with parameters \code{w.penalty}. The default is a
    vector of ones, which reduces to a uniform prior on \code{w}.}

  \item{drop.threshold}{Posterior probability threshold for dropping
    mixture components. Should be a positive number close to zero. If,
    at any point during the optimization, all posterior mixture
    assignment probabilities for a given mixture component \code{k} are
    less than \code{drop.threshold}, the mixture weight for component
    \code{k} is automatically set to zero. Set \code{drop.threshold} to
    zero to disable this behaviour. Setting larger values for
    \code{drop.threshold} may improve computation speed at a small cost
    to numerical accuracy of the final results.}
  
  \item{tol}{Convergence tolerance for co-ordinate ascent updates.}

  \item{maxiter}{Maximum number of co-ordinate ascent iterations.}

  \item{update.order}{Order of the co-ordinate ascent updates for
    fitting the variational approximation. The default is
    \code{update.order = 1:p}, where \code{p} is the number of variables
    (the number of columns of \code{X}).}
  
  \item{verbose}{If \code{verbose = TRUE}, print progress of algorithm
    to console.}
}
\details{
  See \url{https://www.overleaf.com/8954189vvpqnwpxhvhq}.
}
\value{
  An object with S3 class \code{c("varbvsmix","list")}.

  \item{n}{Number of data samples used to fit model.}

  \item{mu.cov}{Posterior mean regression coefficients for covariates,
    including intercept.}

  \item{update.sigma}{If \code{TRUE}, residual variance parameter
    \code{sigma} was fit to data.}

  \item{update.sa}{If \code{TRUE}, mixture variances were fit to data.}

  \item{update.w}{If \code{TRUE}, mixture weights were fit to data.}

  \item{w.penalty}{Penalty used for updating mixture weights.}

  \item{drop.threshold}{Posterior probabiltiy threshold used in the
    optimization procedure for setting mixture weights to zero.}
  
  \item{sigma}{Fitted or user-specified residual variance parameter.}
  
  \item{sa}{User-specified mixture variances.}

  \item{w}{Fitted or user-specified mixture weights.}

  \item{alpha}{Variational estimates of posterior mixture assignent
    probabilities.}

  \item{mu}{Variational estimates of posterior mean coefficients.}

  \item{s}{Variational estimates of posterior variances.}

  \item{lfsr}{Local false sign rate (LFSR) for each variable computed
    from variational estimates of posterior assignment probabilities and
    posterior means and variances. See Stephens (2017) for a definition
    of the LFSR.}
  
  \item{logZ}{Variational lower bound to marginal log-likelihood at each
    iteration of the co-ordinate ascent algorithm.}

  \item{err}{Maximum difference in the variational posterior
    probabilities at each iteration of the co-ordinate ascent
    algorithm.}

  \item{nzw}{Number of nonzero mixture components (including the
    "spike") at each iteration of the co-ordinate ascent algorithm.}
}
\references{
  M. Stephens (2017). False discovery rates: a new deal.
  \emph{Biostatistics} \bold{18}, 275--294.
}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{\link{varbvs}}}
\examples{

# Generate the data set.
set.seed(1)
n    <- 200
p    <- 500
X    <- randn(n,p)
sd   <- c(0,0.2,0.5)
w    <- c(0.9,0.05,0.05)
k    <- sample(length(w),p,replace = TRUE,prob = w)
beta <- sd[k] * rnorm(p)
y    <- c(X \%*\% beta + rnorm(n))

# Fit the model to the data, in which the variances of the mixture
# prior are automatically selected.
fit1 <- varbvsmix(X,NULL,y)

# Fit the model, but use only 3 mixture components in the prior
# instead of the default of 20.
fit2 <- varbvsmix(X,NULL,y,3)

# Use the "ground-truth" prior variances (the ones used to simulate
# the data).
fit3 <- varbvsmix(X,NULL,y,sd^2)

# Compare predicted outcomes against observed outcomes.
y.fit1 <- predict(fit1,X)
print(cor(y,y.fit1))

\dontrun{
library(lattice)
print(xyplot(beta.est ~ beta.true,
             data.frame(beta.true = beta,
                        beta.fitted = rowSums(fit$alpha * fit$mu)),
             pch = 20,col = "royalblue",cex = 1))
}
}
