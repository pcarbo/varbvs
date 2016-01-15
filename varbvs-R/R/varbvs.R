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

  # (1) CHECK INPUTS
  # ----------------
  # If input Z is not NULL, it must have as many rows as X.
  if (!is.null(Z) & nrow(Z) != n)
    stop("Inputs X and Z do not match")

  # Add intercept.
  Z <- cbind(1,Z)
  colnames(Z)[1] <- "intercept"
  
  # Input y must have n entries.
  if (length(y) != n)
    stop("Inputs X and y do not match")
  y <- c(y)
  
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
    if (family == "binomial")
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
  if (!is.matrix(logodds)) {
    prior.same <- TRUE
    logodds    <- t(matrix(logodds))
  } else if (nrow(logodds) != p) {
    prior.same <- TRUE
    logodds    <- t(matrix(logodds))
  } else
    prior.same <- FALSE

  # Here is where I ensure that the numbers of candidate hyperparameter
  # settings are consistent.
  ns <- max(length(sigma),length(sa),ncol(logodds))
  if (length(sigma) == 1)
    sigma <- rep(sigma,ns)
  if (length(sa) == 1)
    sa <- rep(sa,ns)
  if (length(sigma) != ns | length(sa) != ns | ncol(logodds) != ns)
    stop("options.sigma, options.sa and options.logodds are inconsistent")

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
    alpha <- rand(p,ns)
    alpha <- alpha / rep.row(colSums(alpha),p)
  } else
    initialize.params.default <- FALSE
  if (nrow(alpha) != p)
    stop("Input alpha must have as many rows as X has columns")
  if (ncol(alpha) == 1)
    alpha <- rep.col(alpha,ns)

  # Set initial estimates of variational parameter mu.
  if (is.null(mu))
    mu <- randn(p,ns)
  else
    initialize.params.default <- FALSE    
  if (nrow(mu) != p)
    error("Input mu must have as many rows as X has columns")
  if (ncol(mu) == 1)
    mu <- rep.col(mu,ns)

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
    eta <- rep.col(eta,ns)

  # Determine whether to update the variational parameter eta. Note this
  # option is only relevant for logistic regression.
  if (is.null(optimize.eta))
    optimize.eta <- optimize.eta.default

  # (3) PREPROCESSING STEPS
  # -----------------------
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior; see Chipman, George and
  # McCulloch, "The Practical Implementation of Bayesian Model
  # Selection," 2001.
  if (family == "gaussian") {
    if (ncol(Z) == 1) {
      X <- X - rep.row(colMeans(X),p)
      y <- y - mean(y)
    } else {

      # The equivalent expressions in MATLAB are  
      #
      #   y = y - Z*((Z'*Z)\(Z'*y))
      #   X = X - Z*((Z'*Z)\(Z'*X))  
      #
      # This should give the same result as centering the columns of X
      # and subtracting the mean from y when we have only one
      # covariate, the intercept.
      y <- y - c(Z %*% solve(crossprod(Z),t(y %*% Z)))
      X <- X - Z %*% solve(crossprod(Z),t(Z) %*% X)
    }
  }

  # Provide a brief summary of the analysis.
  if (verbose) {
    cat("Welcome to           ")
    cat("--       *                              *               \n")
    cat("VARBVS version 0.9.0 ")
    cat("--       |              |               |               \n")
    cat("large-scale Bayesian ")
    cat("--       ||           | |    |          || |     |   |  \n")
    cat("variable selection   ")
    cat("-- |     || | |    |  | ||  ||        |||| ||    |   || \n")
    cat("*********************")
    cat("*******************************************************\n")
    cat("Fitting variational approximation for Bayesian variable",
        "selection model.\n")
    cat(sprintf("family:     %-8s",family))
    cat("   num. hyperparameter settings:",length(sa),"\n")
    cat(sprintf("samples:    %-6d",n))
    cat(sprintf("     convergence tolerance         %0.1e\n",tol))
    cat(sprintf("variables:  %-6d",p))
    cat("     iid variable selection prior:",tf2yn(prior.same),"\n")
    cat(sprintf("covariates: %-6d",max(0,ncol(Z) - 1)))
    cat("     fit prior var. of coefs (sa):",tf2yn(update.sa),"\n")
    cat("intercept:  yes        ")
    if (family == "gaussian")
      cat("fit residual var. (sigma):   ",tf2yn(update.sigma),"\n")
    else if (family == "binomial")
      cat("fit approx. factors (eta):   ",tf2yn(optimize.eta),"\n")
  }

  # (4) INITIALIZE STORAGE FOR THE OUTPUTS
  # --------------------------------------
  # Initialize storage for the variational estimate of the marginal
  # log-likelihood for each hyperparameter setting (logw), and the
  # variances of the regression coefficients (s).
  logw <- rep(0,ns)
  s    <- matrix(0,p,ns)
    
  # (5) FIT BAYESIAN VARIABLE SELECTION MODEL TO DATA
  # -------------------------------------------------
  if (ns == 1) {
    
    # Find a set of parameters that locally minimize the K-L
    # divergence between the approximating distribution and the exact
    # posterior.
    if (verbose) {
      cat("        variational    max.   incl variance params\n")
      cat(" iter   lower bound  change   vars   sigma      sa\n")
    }
    #
    # TO DO: FIX THIS.
    #
    # [logw sigma sa alpha mu s eta] = ...
    #     outerloop(X,Z,y,family,sigma,sa,logodds,alpha,mu,eta,tol,maxiter,...
    #               verbose,[],update_sigma,update_sa,optimize_eta,n0,sa0);
    if (verbose)
      cat("\n")
  } else {

    # If a good initialization isn't already provided, find a good
    # initialization for the variational parameters. Repeat for each
    # candidate setting of the hyperparameters.
    if (initialize.params) {
      if (verbose) {
        cat("Finding best initialization for",ns,"combinations of",
            "hyperparameters.\n");
        cat("-iteration-   variational    max.   incl variance params\n");
        cat("outer inner   lower bound  change   vars   sigma      sa\n");
      }

      # Repeat for each setting of the hyperparameters.
      for (i in 1:ns) {
        # [logw(i) sigma(i) sa(i) alpha(:,i) mu(:,i) s(:,i) eta(:,i)] = ...
        #     outerloop(X,Z,y,family,sigma(i),sa(i),logodds(:,i),alpha(:,i),...
        #               mu(:,i),eta(:,i),tol,maxiter,verbose,i,update_sigma,...
        #               update_sa,optimize_eta,n0,sa0)
      }
      if (verbose)
        cat("\n")

      # Choose an initialization common to all the runs of the
      # coordinate ascent algorithm. This is chosen from the
      # hyperparameters with the highest variational estimate of the
      # marginal likelihood.
      i     <- which.max(logw)
      alpha <- rep.col(alpha[,i],ns)
      mu    <- rep.col(mu[,i],ns)
      eta   <- rep.col(eta[,i],ns)
    }

    # Compute a variational approximation to the posterior distribution
    # for each candidate setting of the hyperparameters.
    cat("Computing marginal likelihood for",ns,"combinations of",
        "hyperparameters.\n")
    cat("-iteration-   variational    max.   incl variance params\n")
    cat("outer inner   lower bound  change   vars   sigma      sa\n")

    # For each setting of the hyperparameters, find a set of
    # parameters that locally minimize the K-L divergence between the
    # approximating distribution and the exact posterior.
    for (i in 1:ns) {
      # [logw(i) sigma(i) sa(i) alpha(:,i) mu(:,i) s(:,i) eta(:,i)] = ...
      #     outerloop(X,Z,y,family,sigma(i),sa(i),logodds(:,i),alpha(:,i),...
      #               mu(:,i),eta(:,i),tol,maxiter,verbose,i,update_sigma,...
      #               update_sa,optimize_eta,n0,sa0)
    }
    if (verbose)
      cat("\n")
  }
    
  # (6) CREATE FINAL OUTPUT
  # ----------------------
  if (family =='gaussian') {
    fit <- list(family = family,ncov = ncol(Z) - 1,n = n,n0 = n0,sa0 = sa0,
                update.sigma = update.sigma,update.sa = update.sa,
                prior.same = prior.same,logw = logw,sigma = sigma,sa = sa,
                logodds = logodds,alpha = alpha,mu = mu,s = s)
  }

  # NOTE: Do some final post-processing here, such as adding names of
  # variables and covariates to data structures.
  if (prior.same)
    fit.logodds <- c(fit.logodds)
  
  return(fit)
}

# ----------------------------------------------------------------------
# This function implements one iteration of the "outer loop".
outerloop <- function (X, Z, y, family, sigma, sa, logodds, alpha, mu, eta,
                       tol, maxiter, verbose, outer.iter, update.sigma,
                       update.sa, optimize.eta, n0, sa0) {
  p <- length(alpha)
  if (length(logodds) == 1)
    logodds <- rep(logodds,p)
  if (family == "gaussian") {
    # TO DO: FIX THIS.
    return(varbvsnorm(X,y,sigma,sa,log(10)*logodds,alpha,mu,tol,maxiter,
                      verbose,outer.iter,update_sigma,update.sa,n0,sa0))
  } else if (family == "binomial" & ncol(Z) == 1) {
    # [logw sa alpha mu s eta] = ...
    #     varbvsbin(X,y,sa,log(10)*logodds,alpha,mu,eta,tol,maxiter,verbose,...
    #               outer_iter,update_sa,optimize_eta,n0,sa0)
  } else if (family == "binomial") {
    # [logw sa alpha mu s eta] = ...
    #     varbvsbinz(X,Z,y,sa,log(10)*logodds,alpha,mu,eta,tol,maxiter,...
    #               verbose,outer_iter,update_sa,optimize_eta,n0,sa0)
  }
  return(out)
}
