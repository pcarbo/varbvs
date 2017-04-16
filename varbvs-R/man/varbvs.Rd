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
  varbvs (X, Z, y, family = c("gaussian","binomial"), sigma, sa,
          logodds, alpha, mu, eta, update.sigma, update.sa,
          optimize.eta, initialize.params, nr = 100, sa0 = 0, n0 = 0,
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
  "binomial"}. If missing, residual variance parameter is automatically
  fitted to data by computing approximate maximum-likelihood (ML)
  estimate.}

  \item{sa}{Hyperparameter sa is the prior variance of regression
    coefficients for variables that are included in the model. This
    prior variance is always scaled by sigma (for logistic regression,
    we take sigma = 1). Scaling the variance of the coefficients in this
    way is necessary to ensure that this prior is invariant to
    measurement scale (e.g., switching from grams to kilograms). This
    input specifies the candidate settings for sa, of the same length as
    inputs sigma and logodds (or have length equal to the number of
    columns of logodds). If missing, prior variance is automatically
    fitted to data by compute approximate maximum (ML) estimates, or
    maximum a posteriori estimates when \code{n0 > 0} and \code{sa0 >
    0}.}

  \item{logodds}{Hyperparameter logodds is the prior log-odds that a
    variable is included in the regression model; it is defined as
    \eqn{logodds = log10(q/(1-q)),} where q is the prior probability
    that a variable is included in the regression model. Note that we
    use the base-10 logarithm instead of the natural logarithm because
    it is usually more natural to specify prior log-odds settings in
    this way. The prior log-odds may also be specified separately for
    each variable, which is useful is there is prior information about
    which variables are most relevant to the outcome. This is
    accomplished by setting logodds to a p x ns matrix, where p is the
    number of variables, and ns is the number of hyperparameter
    settings. Note it is not possible to fit the logodds parameter; if
    logodds input is not provided as input, then it is set to the
    default value when sa and sigma are missing, and otherwise an error
    is generated.}

  \item{alpha}{Good initial estimate for the variational parameter alpha
    for each hyperparameter setting. Either missing, or a p x ns matrix,
    where p is the number of variables, and ns is the number of
    hyperparameter settings.}

  \item{mu}{Good initial estimate for the variational parameter mu for
    each hyperparameter setting. Either missing, or a p x ns matrix,
    where p is the number of variables, and ns is the number of
    hyperparameter settings.}

  \item{eta}{Good initial estimate of the additional free parameters
    specifying the variational approximation to the logistic regression
    factors. Either missing, or an n x ns matrix, where n is the number
    of samples, and ns is the number of hyperparameter settings.}

  \item{update.sigma}{Setting this to TRUE ensures that sigma is always
    fitted to data, in which case input vector sigma is used to provide
    initial estimates.}

  \item{update.sa}{Setting this to TRUE ensures that sa is always fitted
    to data, in which case input vector sa is used to provide initial
    estimates.}

  \item{optimize.eta}{When optimize.eta = TRUE, eta is fitted to the
    data during the inner loop coordinate ascent updates, even when good
    estimates of eta are provided as input.}

  \item{initialize.params}{If TRUE, the initialization stage of the
    variational inference algorithm (see below) will be skipped, which
    saves computation time for large data sets.}

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
\section{Regression models}{
    Two types of outcomes (y) are modeled: (1) a continuous outcome,
    also a "quantitative trait" in the genetics literature; or (2) a
    binary outcome with possible values 0 and 1. For the former, set
    \code{family = "gaussian"}, in which case, the outcome is
    i.i.d. normal with mean \eqn{u0 + Z*u + X*b} and variance sigma, in
    which u and b are vectors of regresion coefficients, and u0 is the
    intercept. In the second case, we use logistic regression to model
    the outcome, in which the probability that y = 1 is equal to
    \eqn{sigmoid(u0 + Z*u + X*b).} See \code{help(sigmoid)} for a
    description of the sigmoid function. Note that the regression always
    includes an intercept term (u0).
}
\section{Co-ordinate ascent optimization procedure}{
    For both regression models, the fitting procedure consists of an inner
    loop and an outer loop. The outer loop iterates over each of the
    hyperparameter settings (sa, sigma and logodds). Given a setting of
    the hyperparameters, the inner loop cycles through coordinate ascent
    updates to tighten the lower bound on the marginal likelihood,
    \eqn{Pr(Y | X, sigma, sa, logodds)}. The inner loop coordinate
    ascent updates terminate when either (1) the maximum number of inner
    loop iterations is reached, as specified by \code{maxiter}, or (2)
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
}
\section{Variational approximation}{
    Outputs alpha, mu and s specify the approximate posterior distribution
    of the regression coefficients. Each of these outputs is a p x ns
    matrix. For the ith hyperparameter setting, alpha[,i] is the
    variational estimate of the posterior inclusion probability (PIP)
    for each variable; mu[,i] is the variational estimate of the
    posterior mean coefficient given that it is included in the model;
    and s[,i] is the estimated posterior variance of the coefficient
    given that it is included in the model. These are also the
    quantities that are optimized as part of the inner loop coordinate
    ascent updates. An additional free parameter, eta, is needed for
    fast computation with the logistic regression model \code{(family =
    "binomial")}. The fitted value of eta is returned as an n x ns
    matrix.

    The variational estimates should be interpreted carefully,
    especially when variables are strongly correlated. For example,
    consider the simple scenario in which 2 candidate variables are
    closely correlated, and at least one of them explains the outcome
    with probability close to 1. Under the correct posterior
    distribution, we would expect that each variable is included with
    probability ~0.5. However, the variational approximation, due to the
    conditional independence assumption, will typically get this wrong,
    and concentrate most of the posterior weight on one variable (the
    actual variable that is chosen will depend on the starting
    conditions of the optimization). Although the individual PIPs are
    incorrect, a statistic summarizing the variable selection for both
    correlated variables (e.g., the total number of variables included
    in the model) should be reasonably accurate.

    More generally, if variables can be reasonably grouped together
    based on their correlations, we recommend interpreting the variable
    selection results at a group level. For example, in genome-wide
    association studies (see the vignettes) ,a SNP with a high PIP
    indicates that this SNP is probably associated with the trait, and
    one or more nearby SNPs within a chromosomal region, or ``locus,''
    may be associated as well. Therefore, we interpreted the GWAS
    variable selection results at the level of loci, rather than at the
    level of individual SNPs. 
    
    Also note that special care is required for interpreting the results
    of the variational approximation with the logistic regression
    model. In particular, interpretation of the individual estimates of
    the regression coefficients (e.g., the posterior mean estimates
    \code{fit$mu}) is not straightforward due to the additional
    approximation introduced on the individual nonlinear factors in the
    likelihood. As a general guideline, only the relative magnitudes of
    the coefficients are meaningful.
}
\section{Averaging over hyperparameter settings}{
    In many settings, it is good practice to account for uncertainty in
    the hyperparameters when reporting final posterior quantities. For
    example, hyperparameter sa is often estimated with a high degree of
    uncertainty when only a few variables are included in the model.
    Provided that (1) the hyperparameter settings sigma, sa and logodds
    adequately represent the space of possible hyperparameter settings
    with high posterior mass, (2) the hyperparameter settings are drawn
    from the same distribution as the prior, and (3) the
    fully-factorized variational approximation closely approximates the
    true posterior distribution, then final posterior quantities can be
    calculated by using logw as (unnormalized) log-marginal
    probabilities.

    Even when conditions (1), (2) and/or (3) are not satisfied, this can
    approach can still often yield reasonable estimates of averaged
    posterior quantities. The examples below demonstrate how final
    posterior quantities are reported by function
    \code{\link{summary.varbvs}} (see \code{help(summary.varbvs)} for
    more details). To account for discrepancies between the prior on
    (sigma,sa,logodds) and the sampling density used to draw candidate
    settings of the hyperparameters, adjust the log-probabilities by
    setting \code{fit$logw <- fit$logw + logp/logq}, where logp is the
    log-density of the prior distribution, and logq is the log-density
    of the sampling distribution. (This is importance sampling; see, for
    example, R. M. Neal, Annealed importance sampling, \emph{Statistics
      and Computing}, 2001.)
}
\section{Prior on proportion of variance explained}{
  Specifying the prior variance of the regression coefficients (sa) can
  be difficult, which is why we have included the option of fitting this
  hyperparameter to the data (see input argument \code{update.sa}
  above). However, in many settings, especially when a small number of
  variables are included in the regression model, it is preferrable to
  average over candidate settings of sa instead of fitting sa to the
  data. To choose a set of candidate settings for sa, we have advocated
  for setting sa indirectly through a prior estimate of the proportion
  of variance in the outcome explained by the variables (abbreviated as
  PVE), since it is often more natural to specify the PVE rather than
  the prior variance (see references below). This is technically only
  suitable or the linear regression model (\code{family = "gaussian"}),
  but could potentially be used for the linear regression model in an
  approximate way.

  For example, one could approximate a uniform prior on the PVE by
  drawing the PVE uniformly between 0 and 1, additionally specifying
  candidate settings for the prior log-odds, then computing the prior
  variance (sa) as follows:\preformatted{sx <- sum(var1.cols(X))
sa <- PVE/(1-PVE)/(sigmoid(log(10)*logodds)*sx)}
  
  Note that this calculation will yield \code{sa = 0} when \code{PVE =
  0}, and \code{sa = Inf} when \code{PVE = 1}.
  
  Also, bear in mind that if there are additional covariates (Z)
  included in the linear regression model that explain variance in Y,
  then it will usually make more sense to first remove the linear
  effects of these covariates before performing this calculation. The
  PVE would then represent the prior proportion of variance in the
  residuals of Y that are explained by the candidate variables. For an
  example of how to do this, see the code for function \code{varbvs}, in
  file varbvs.R, under "preprocessing steps". Alternatively, one could
  include the matrix Z in the calculation above, taking care to ensure
  that the covariates are included in the model with probability 1.
}
\section{Memory requirements}{
    Finally, we point out that the optimization procedures were
    carefully designed so that they can be applied to very large data
    sets; to date, this code has been tested on data sets with >500,000
    variables and >10,000 samples. An important limiting factor is the
    ability to store the data matrix X in memory. To reduce memory
    requirements, in the MATLAB interface we require that X be single
    precision, but this option is not available in R. Additionally, we
    mostly avoid generating intermediate products that are of the same
    size as X. Only one such intermediate product is generated when
    \code{family = "gaussian"}, and none for \code{family = "binomial"}.
}
\value{
    An object with S3 class \code{c("varbvs","list")}.

  \item{family}{Either "gaussian" or "binomial".}

  \item{n}{Number of data samples used to fit model.}

  \item{ncov}{Number of covariates (columns of Z).}

  \item{sigma}{Settings for sigma (family = "gaussian" only).}

  \item{sa}{Settings for prior variance parameter.}

  \item{logodds}{Prior log-odds settings.}

  \item{prior.same}{TRUE if prior is identical for all variables. When
        logodds is a p x ns matrix, \code{prior.same = FALSE}.}

  \item{sa0}{Scale parameter for prior on hyperparameter sa.}

  \item{n0}{Degrees of freedom for prior on hyperparameter sa.}

  \item{update.sigma}{If TRUE, sigma was fit to data for each setting of
    prior logodds (\code{family = "gaussian"} only).}

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

  \item{mu.cov}{Posterior mean regression coefficients for covariates,
  including intercept, for each hyperparameter setting.}
  
  \item{eta}{Additional variational parameters for \code{family =
    "binomial"} only.}
  
  \item{optimize.eta}{If TRUE, eta was fit to data (\code{family =
  "binomial"} only).}

  \item{pve}{For each hyperparameter setting, and for each variable,
    mean estimate of the proportion of variance in outcome explained
    conditioned on variable being included in the model (\code{family =
    "gaussian"} only).}
    
  \item{model.pve}{Samples drawn from posterior distribution giving
  estimates of proportion of variance in outcome (y) explained by fitted
  variable selection model. This is for \code{family = "gaussian"}
  only). To obtain the posterior mean estimate of the proportion of
  variance explained (PVE), for example, simply type
  \code{mean(fit$model.pve)}.}
}
\references{
  P. Carbonetto and M. Stephens (2012). Scalable variational inference
  for Bayesian variable selection in regression, and its accuracy in
  genetic association studies. \emph{Bayesian Analysis} \bold{7},
  73--108.

  Y. Guan and M. Stephens (2011). Bayesian variable selection regression
  for genome-wide association studies and other large-scale
  problems. \emph{Annals of Applied Statistics} \bold{5}, 1780--1815.

  X. Zhou, P. Carbonetto and M. Stephens (2013). Polygenic modeling with
  Bayesian sparse linear mixed models. \emph{PLoS Genetics} \bold{9},
  e1003264.
}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{summary.varbvs}, \code{varbvscoefcred}, \code{varbvspve},
  \code{varbvsnorm}, \code{varbvsbin}, \code{varbvsbinz},
  \code{normalizelogweights}, \code{varbvs-internal}
}
\examples{

# LINEAR REGRESSION EXAMPLE
# -------------------------
# Data are 200 uncorrelated ("unlinked") single nucleotide polymorphisms
# (SNPs) with simulated genotypes, in which the first 20 of them have an
# effect on the outcome. Also generate data for 3 covariates.
maf <- 0.05 + 0.45*runif(200)
X   <- (runif(400*200) < maf) + (runif(400*200) < maf)
X   <- matrix(as.double(X),400,200,byrow = TRUE)
Z   <- randn(400,3)

# Generate the ground-truth regression coefficients for the variables
# (X) and additional 3 covariates (Z). Adjust the QTL effects so that
# the variables (SNPs) explain 50 percent of the variance in the
# outcome.
u    <- c(-1,2,1)
beta <- c(rnorm(20),rep(0,180))
beta <- 1/sd(c(X \%*\% beta)) * beta

# Generate the quantitative trait measurements.
y <- c(-2 + Z \%*\% u + X \%*\% beta + rnorm(400))

# Fit the variable selection model.
fit <- varbvs(X,Z,y,logodds = seq(-3,-1,0.1))
print(summary(fit))

# Compute posterior inclusion probabilities (pip) and posterior mean
# regression coefficients (beta) averaged over the hyperparameter
# settings.
w    <- normalizelogweights(fit$logw);
pip  <- fit$alpha \%*\% c(w)
beta <- fit$mu \%*\% c(w)

# Compute the posterior mean estimate of hyperparameter sa.
sa <- sum(fit$sa * w)

# Compare estimated outcomes against observed outcomes.
y.fit <- predict(fit,X,Z)
print(cor(y,y.fit))

# LOGISTIC REGRESSION EXAMPLE
# ---------------------------
# Data are 100 uncorrelated ("unlinked") single nucleotide polymorphisms
# (SNPs) with simulated genotypes, in which the first 10 of them have an
# effect on the outcome. Also generate data for 2 covariates.
maf <- 0.05 + 0.45*runif(100)
X   <- (runif(750*100) < maf) + (runif(750*100) < maf)
X   <- matrix(as.double(X),750,100,byrow = TRUE)
Z   <- randn(750,2)

# Generate the ground-truth regression coefficients for the variables
# (X) and additional 2 covariates (Z).
u    <- c(-1,1)
beta <- c(0.5*rnorm(10),rep(0,90))

# Simulate the binary trait (case-control status) as a coin toss with
# success rates given by the logistic regression.
y <- as.double(runif(750) < sigmoid(-1 + Z \%*\% u + X \%*\% beta))

# Fit the variable selection model.
fit <- varbvs(X,Z,y,"binomial",logodds = seq(-2,-0.5,0.5))
print(summary(fit))
}
