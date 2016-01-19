\name{varbvs}
\alias{varbvs}
\title{Fit variable selection model using variational approximation methods.}

\description{Compute fully-factorized variational approximation for
  Bayesian variable selection in linear (family = gaussian) or logistic
  regression (family = binomial). More precisely, find the "best"
  fully-factorized approximation to the posterior distribution of the
  coefficients, with spike-and-slab priors on the coefficients. By
  "best", we mean the approximating distribution that locally minimizes
  the Kullback-Leibler divergence between the approximating distribution
  and the exact posterior.
}
\usage{
  varbvs(X, Z, y, family = "gaussian", sigma = NULL, sa = NULL,
         logodds = NULL, alpha = NULL, mu = NULL, eta = NULL,
         update.sigma = NULL, update.sa = NULL, optimize.eta = NULL,
         initialize.params = NULL, nr = 100, sa0 = 0, n0 = 0,
         tol = 1e-4, maxiter = 1e4, verbose = TRUE)
}
\arguments{
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

  \item{family}{"gaussian" for linear regression model, or "binomial"
                for logistic regression model.}

  \item{sigma}{Candidate settings for the residual variance
  parameter. Must be of the same length as inputs sa and logodds (or
  have length equal to the number of columns of logodds). Only used for
  linear regression, and will generate an error if \code{family =
  "binomial"}. If sigma = NULL, residual variance parameter is
  automatically fitted to data by computing approximate
  maximum-likelihood (ML) estimate.}

  \item{sa}{Hyperparameter sa is the prior variance of regression
    coefficients for variables that are included in the model. This
    input specifies the candidate settings for sa, of the same length as
    inputs sigma and logodds (or have length equal to the number of
    columns of logodds). If sa = NULL, prior variance is automatically
    fitted to data by compute approximate maximum (ML) estimates, or
    maximum a posteriori estimates when \code{n0 > 0} and \code{sa0 > 0}.}

  \item{logodds}{Hyperparameter logodds is the prior log-odds that a
    variable is included in the regression model; it is defined as
    \deqn{logodds = log10(q/(1-q)),} where q is the prior probability
    that a variable is included in the regression model. The prior
    log-odds may also be specified separately for each variable, which
    is useful is there is prior information about which variables are
    most relevant to the outcome. This is accomplished by setting
    logodds to a p x ns matrix, where p is the number of variables, and
    ns is the number of hyperparameter settings. In this case,
    \code{prior.same = false} in the return value. Note it is not
    possible to fit the logodds parameter; if logodds input is not
    provided as input, then it is set to the default value when sa and
    sigma are NULL, and otherwise an error is generated.}

  \item{alpha}{Good initial estimate for the variational parameter
  alpha for each hyperparameter setting. Either NULL, or a p x ns
  matrix, where p is the number of variables, and ns is the number of
  hyperparameter settings.}

  \item{mu}{Good initial estimate for the variational parameter mu for
  each hyperparameter setting. Either NULL, or a p x ns matrix, where p
  is the number of variables, and ns is the number of hyperparameter
  settings.}

  \item{eta}{Good initial estimate of the additional free parameters
  specifying the variational approximation to the logistic regression
  factors. Either NULL, or an n x ns matrix, where n is the number of
  samples, and ns is the number of hyperparameter settings.}

  \item{update.sigma}{Setting this to TRUE ensures that sigma is always
  fitted to data, in which case input vector sigma is used to provide
  initial estimates.}

  \item{update.sa}{Setting this to TRUE ensures that sa is always
  fitted to data, in which case input vector sa is used to provide
  initial estimates.}

  \item{optimize.eta}{When optimize.eta = TRUE, eta is fitted to the
  data during the inner loop coordinate ascent updates.}

  \item{initialize.params}{If TRUE, the initialization stage of the
    variational inference algorithm (see below) will be skipped.}

  \item{nr}{Number of samples of "model PVE" to draw from posterior.}

  \item{sa0}{Scale parameter for a scaled inverse chi-square prior on
  hyperparameter sa. Must be >= 0.}
  
  \item{n0}{Number of degrees of freedom for a scaled inverse chi-square
  prior on hyperparameter sa. Must be >= 0.}

  \item{tol}{Convergence tolerance for inner loop.}

  \item{maxiter}{Maximum number of inner loop iterations.}

  \item{verbose}{If verbose = TRUE, print progress of algorithm to
    console.}
  
}
\details{
  
    \bold{Regression models.} Two types of outcomes (y) are modeled: (1)
    a continuous outcome, also a "quantitative trait" in the genetics
    literature; or (2) a binary outcome with possible values 0 and
    1. For the former, set \code{family = "gaussian"}, in which case,
    the outcome is i.i.d. normal with mean \deqn{u0 + Z*u + X*b} and
    variance sigma, in which u and b are vectors of regresion
    coefficients, and u0 is the intercept. In the second case, we use
    logistic regression to model the outcome, in which the probability
    that y = 1 is equal to \deqn{sigmoid(u0 + Z*u + X*b).} See
    \code{help(sigmoid)} for a description of the sigmoid function. Note
    that the regression always includes an intercept term (u0).

    \bold{Co-ordinate ascent optimization procedure.} For both
    regression models, the fitting procedure consists of an inner loop
    and an outer loop. The outer loop iterates over each of the
    hyperparameter settings (sa, sigma and logodds). Given a setting of
    the hyperparameters, the inner loop cycles through coordinate ascent
    updates to tighten the lower bound on the marginal likelihood,
    \deqn{Pr(Y | X, sigma, sa, logodds)}. The inner loop coordinate
    ascent updates terminate when either (1) the maximum number of inner
    loop iterations is reached, as specified by options.maxiter, or (2)
    the maximum difference between the estimated posterior inclusion
    probabilities is less than \code{tol}.

    To provide a more accurate variational approximation of the posterior
    distribution, by default the fitting procedure has two stages. In the
    first stage, the entire fitting procedure is run to completion, and the
    variational parameters (alpha, mu, s, eta) corresponding to the maximum
    lower bound are then used to initialize the coordinate ascent updates
    in a second stage. Although this has the effect of doubling the
    computation time (in the worst case), the final posterior estimates
    tend to be more accurate with this two-stage fitting procedure. 
    
    \bold{Variational approximation.} Outputs alpha, mu and s specify
    the approximate posterior distribution of the regression
    coefficients. Each of these outputs is a p x ns matrix. For the ith
    hyperparameter setting, alpha[,i] is the variational estimate of the
    posterior inclusion probability (PIP) for each variable; mu[,i] is
    the variational estimate of the posterior mean coefficient given
    that it is included in the model; and s[,i] is the estimated
    posterior variance of the coefficient given that it is included in
    the model. These are also the quantities that are optimized as part
    of the inner loop coordinate ascent updates. An additional free
    parameter, eta, is needed for fast computation with the logistic
    regression model \code{(family = "binomial")}. The fitted value of
    eta is returned as an n x ns matrix.

    \bold{Correctly averaging over hyperparameter settings.} In many
    settings, it is good practice to account for uncertainty in the
    hyperparameters when reporting final posterior quantities. For
    example, hyperparameter sa is often estimated with a high degree of
    uncertainty when only a few variables are included in the model.
    Provided that (1) the hyperparameter settings sigma, sa and logodds
    adequately represent the space of possible hyperparameter settings
    with high posterior mass, (2) the hyperparameter settings are drawn
    from the same distribution as the prior, and (3) the
    fully-factorized variational approximation closely approximates the
    true posterior distribution, then final posterior quantities can be
    calculated by using logw as (unnormalized) log-importance weights.
    Even when conditions (1), (2) and/or (3) are not satisfied, this can
    approach can still often yield reasonable estimates of averaged
    posterior quantities. For example, do the following to compute
    posterior inclusion probabilities (PIPs) averaged over the
    hyperparameter settings:
    \code{
      w <- normalizelogweights(fit$logw);\cr
      PIP <- fit$alpha %*% c(w)
    }
    And, for example, do the following to compute the posterior mean
    estimate of hyperparameter sa:
    \code{
        w <- normalizelogweights(fit$logw)
        mean_sa <- dot(fit$sa,w)
    }
    This is precisely how final posterior quantities are reported by
    function \code{varbvsprint} (see \code{help(varbvsprint)}for more
    details). To account for discrepancies between the prior on
    (sigma,sa,logodds) and the sampling density used to draw candidate
    settings of the hyperparameters, adjust the log-importance weights
    by setting \code{fit$logw <- fit$logw + logp/logq}, where logp is the
    log-density of the prior distribution, and logq is the log-density
    of the sampling distribution.

    \bold{Memory requirements.} Finally, we point out that the
    optimization procedures were carefully designed so that they can be
    applied to very large data sets; to date, this code has been tested
    on data sets with >500,000 variables and >10,000 samples. An
    important limiting factor is the ability to store the data matrix X
    in memory. To reduce memory requirements, in the MATLAB interface we
    require that X be single precision, but this option is not available
    in R. Additionally, we mostly avoid generating intermediate products
    that are of the same size as X. Only one such intermediate product
    is generated when \code{family = "gaussian"}, and none for
    \code{family = "binomial"}.
    
}  
\value{

  An object with S3 class \code{"varbvs","list"}.

  \item{family}{Either "gaussian" or "binomial".}

  \item{n}{Number of data samples used to fit model.}

  \item{ncov}{Number of covariates (columns of Z).}

  \item{sigma}{Settings for sigma (family = "gaussian" only).}

  \item{sa}{Settings for prior variance parameter.}

  \item{logodds}{Prior log-odds settings.}

  \item{prior.same}{TRUE if prior is identical for all variables.}

  \item{sa0}{Scale parameter for prior on hyperparameter sa.}

  \item{n0}{Degrees of freedom for prior on hyperparameter sa.}

  \item{update.sigma}{If TRUE, sigma was fit to data for each setting of
    prior logodds (family = "gaussian" only).}

  \item{update.sa}{If TRUE, sa was fit to data for each setting of prior
  logodds.}

  \item{logw}{An array with ns elements, in which \code{logw[i]} is the
    variational lower bound on the marginal log-likelihood for setting i
    of the hyperparameters.}

  \item{alpha}{Variational estimates of posterior inclusion
    probabilities for each hyperparameter setting.}

  \item{mu}{Variational estimates of posterior mean coefficients for
    each hyperparameter setting.}

  \item{s}{Variational estimates of posterior variances for each
    hyperparameter setting.}

  \item{eta}{Additional variational parameters for family = "binomial"
    only.}
  
  \item{eta}{If TRUE, eta was fit to data (family = 'binomial' only).}

  \item{pve}{For each hyperparameter setting, and for each variable,
    mean estimate of the proportion of variance in outcome explained
    conditioned on variable being included in the model (family =
    "gaussian" only).}

  \item{model.pve}{Samples drawn fom posterior of proportion of variance
  in outcome explained by variable selection model (only for family =
  "gaussian").}
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
# Gaussian
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit1=glmnet(x,y)
print(fit1)
coef(fit1,s=0.01) # extract coefficients at a single value of lambda
predict(fit1,newx=x[1:10,],s=c(0.01,0.005)) # make predictions
}

 
