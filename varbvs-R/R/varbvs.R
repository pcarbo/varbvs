# Compute fully-factorized variational approximation for Bayesian
# variable selection in linear (family = "gaussian") or logistic
# regression (family = "binomial"). See varbvs.Rd for details.
varbvs <- function (X, Z, y, family = "gaussian", sigma = NULL, sa = NULL,
                    logodds = NULL, alpha = NULL, mu = NULL, eta = NULL,
                    update.sigma = NULL, update.sa = NULL, optimize.eta = NULL,
                    initialize.params = NULL, sa0 = 0, n0 = 0, tol = 1e-4,
                    maxiter = 1e4, verbose = TRUE) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  browser()
  
  # (1) CHECK INPUTS
  # ----------------
  # If input Z is not NULL, it must have as many rows as X.
  if (!is.null(Z) & nrow(Z) != n)
      stop("Inputs X and Z do not match")

  # Add intercept.
  Z <- cbind(matrix(1,n,1),Z)

  # Input y must have n entries.
  if (length(y) != n)
    stop("Inputs X and y do not match")

  # Check choice of regression model.
  if (family != "gaussian" & family != "binomial")
    stop("family must be gaussian or binomial")

  # (2) PROCESS OPTIONS
  # -------------------
  # Get candidate settings for the variance of the residual (sigma),
  # if provided. Note that this option is not valid for a binary trait.
  if (is.null(sigma)) {
    sigma <- var(y)
    update.sigma.default <- TRUE
  } else {
    sigma <- c(sigma)
    update.sigma.default <- FALSE
    else if (family == "binomial")
      stop("Input sigma is not allowed for family = binomial")
  }
  
  # Get candidate settings for the prior variance of the coefficients
  # (sa), if provided.
  if (is.null(sa)) {
    sa <- 1
    update.sa.default <- TRUE
  } else {
    sa <- c(sa)
    update.sa.default <- FALSE
  }

  # Get candidate settings for the prior log-odds of inclusion. This
  # may either be specified as a vector, in which case this is the
  # prior applied uniformly to all variables, or it is a p x ns
  # matrix, where p is the number of variables and ns is the number of
  # candidate hyperparameter settings, in which case the prior
  # log-odds is specified separately for each variable. A default
  # setting is only available if the number of other hyperparameter
  # settings is 1, in which case we select 20 candidate settings for
  # the prior log-odds, evenly spaced between log10(1/p) and
  # log10(0.5).
  if (is.null(logodds)) {
    if (length(sigma) == 1 & length(sa) == 1)
      logodds <- seq(-log10(p),-0.3,length.out = 20)
    else
      stop("logodds can only be NULL when length(sigma) = length(sa) = 1")
  }
  if (is.matrix(logodds) & nrow(logodds) == p)
    prior.same <- FALSE
  else {
    prior.same <- TRUE
    logodds    <- t(matrix(logodds))
  }

  # Here is where I ensure that the numbers of candidate hyperparameter
  # settings are consistent.
  ns <- max(length(sigma),length(sa),ncol(logodds))
  if (length(sigma) == 1)
    sigma <- rep(sigma,ns)
  if (length(sa) == 1)
    sa <- rep(sa,ns)
  if (length(sigma) != ns | length(sa) != ns | ncol(logodds) != ns)
    stop("options.sigma, options.sa and options.logodds are inconsistent")
  if (nrow(logodds) == 1)
    logodds <- c(logodds)

  # Determine whether to update the residual variance parameter. Note
  # that this option is only relevant for a binary trait.
  if (is.null(update.sigma))
    update.sigma <- update.sigma.default

  # Determine whether to update the prior variance of the regression
  # coefficients.
  if (is.null(update.sa))
    update.sa <- update.sa.default
  
  # Set initial estimates of variational parameter alpha.
  initialize.params.default <- TRUE
  if (is.null(alpha)) {
    alpha <- matrix(runif(p*ns),p,ns)
    alpha <- alpha / colSums(alpha)
  } else
    initialize.params.default <- FALSE
  if (nrow(alpha) != p)
    stop("Input alpha must have as many rows as X has columns")
  if (ncol(alpha) == 1)
    alpha <- matrix(rep(alpha,1,ns),p,ns)

  # Set initial estimates of variational parameter mu.
  if (is.null(mu))
    mu <- matrix(rnorm(p*ns),p,ns)
  else
    initialize.params.default <- FALSE    
  if (nrow(mu) != p)
    error("'Input mu must have as many rows as X has columns")
  if (ncol(mu) == 1)
    mu <- matrix(rep(mu,1,ns),p,ns)

  # Determine whether to find a good initialization for the
  # variational parameters.
  if (is.null(initialize.params))
    initialize.params <- initialize.params.default

  # Set initial estimates of variational parameter eta. Note this
  # input is only relevant for logistic regression.
  if (is.null(eta)) {
    eta                  <- matrix(1,n,ns)
    optimize.eta.default <- TRUE
  } else {
    optimize.eta.default <- FALSE
    if (family != "binomial")
      stop("Input eta is only valid for family = binomial")
  }
  if (nrow(eta) != n)
    stop("Input eta must have as many rows as X")
  if (ncol(eta) == 1)
    eta <- matrix(rep(eta,ns),n,ns)

  # Determine whether to update the variational parameter eta. Note this
  # option is only relevant for logistic regression.
  if (is.null(optimize.eta))
    optimize.eta <- optimize.eta.default
}
