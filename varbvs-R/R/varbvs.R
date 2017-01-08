# Compute fully-factorized variational approximation for Bayesian
# variable selection in linear (family = "gaussian") or logistic
# regression (family = "binomial"). See varbvs.Rd for details.
varbvs <- function (X, Z, y, family = c("gaussian","binomial"), sigma, sa,
                    logodds, alpha, mu, eta, update.sigma, update.sa,
                    optimize.eta, initialize.params, nr = 100, sa0 = 0,
                    n0 = 0, tol = 1e-4, maxiter = 1e4, verbose = TRUE) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # (1) CHECK INPUTS
  # ----------------
  # Check input X.
  if (!(is.matrix(X) & is.double(X) & sum(is.na(X)) == 0))
    stop("Input X must be a double-precision matrix with no missing values.")
  
  # Add column names to X if they are not provided.
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X",1:p)
  
  # Check input Z.
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) | sum(is.na(Z)) > 0)
      stop("Input Z must be a numeric matrix with no missing values.")
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    storage.mode(Z) <- "double"
  }
  
  # Add intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)

  # Check input y.
  if (!is.numeric(y) | sum(is.na(y)) > 0)
    stop("Input y must be a numeric vector with no missing values.")
  if (length(y) != n)
    stop("Inputs X and y do not match")
  y <- c(as.double(y))
  
  # Get choice of regression model.
  family <- match.arg(family)

  # (2) PROCESS OPTIONS
  # -------------------
  if (!is.finite(maxiter))
    stop("Input maxiter must be a finite number")
  
  # Get candidate settings for the variance of the residual (sigma),
  # if provided. Note that this option is not valid for a binary trait.
  if (missing(sigma)) {
    sigma <- var(y)
    update.sigma.default <- TRUE
  } else {
    sigma <- c(sigma)
    update.sigma.default <- FALSE
    if (family == "binomial")
      stop("Input sigma is not allowed for family = binomial")
  }
  
  # Get candidate settings for the prior variance of the coefficients,
  # if provided.
  if (missing(sa)) {
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
  # the prior log-odds, evenly spaced between log10(1/p) and -1.
  if (missing(logodds)) {
    if (length(sigma) == 1 & length(sa) == 1)
      logodds <- seq(-log10(p),-1,length.out = 20)
    else
      stop("logodds can only be missing when length(sigma) = length(sa) = 1")
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
  if (ncol(logodds) == 1)
    logodds <- rep.col(logodds,ns)
  if (length(sigma) != ns | length(sa) != ns | ncol(logodds) != ns)
    stop("options.sigma, options.sa and options.logodds are inconsistent")

  # Determine whether to update the residual variance parameter. Note
  # that this option is only relevant for a binary trait.
  if (missing(update.sigma))
    update.sigma <- update.sigma.default

  # Determine whether to update the prior variance of the regression
  # coefficients.
  if (missing(update.sa))
    update.sa <- update.sa.default
  
  # Set initial estimates of variational parameter alpha.
  initialize.params.default <- TRUE
  if (missing(alpha)) {
    alpha <- rand(p,ns)
    alpha <- alpha / rep.row(colSums(alpha),p)
  } else
    initialize.params.default <- FALSE
  if (nrow(alpha) != p)
    stop("Input alpha must have as many rows as X has columns")
  if (ncol(alpha) == 1)
    alpha <- rep.col(alpha,ns)

  # Set initial estimates of variational parameter mu.
  if (missing(mu))
    mu <- randn(p,ns)
  else
    initialize.params.default <- FALSE    
  if (nrow(mu) != p)
    stop("Input mu must have as many rows as X has columns")
  if (ncol(mu) == 1)
    mu <- rep.col(mu,ns)

  # Determine whether to find a good initialization for the
  # variational parameters.
  if (missing(initialize.params))
    initialize.params <- initialize.params.default
  else if (initialize.params & ns == 1)
    stop(paste("initialize.params = TRUE has no effect when there is",
               "only one hyperparameter setting"))

  # Set initial estimates of variational parameter eta. Note this
  # input is only relevant for logistic regression.
  if (missing(eta)) {
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
  if (missing(optimize.eta))
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

    # Here I compute two quantities that are used here to remove
    # linear effects of the covariates (Z) on X and y, and later on
    # (in function "outerloop"), to efficiently compute estimates of
    # the regression coefficients for the covariates.
    SZy <- solve(crossprod(Z),c(y %*% Z))
    SZX <- solve(crossprod(Z),t(Z) %*% X)
    if (ncol(Z) == 1) {
      X <- X - rep.row(colMeans(X),n)
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
      y <- y - c(Z %*% SZy)
      X <- X - Z %*% SZX
    }
  } else {
    SZy <- NULL
    SZX <- NULL
  }

  # Provide a brief summary of the analysis.
  if (verbose) {
    cat("Welcome to           ")
    cat("--       *                              *               \n")
    cat("VARBVS version 2.0.3 ")
    cat("--       |              |               |               \n")
    cat("large-scale Bayesian ")
    cat("--       ||           | |    |          || |     |   |  \n")
    cat("variable selection   ")
    cat("-- |     || | |    |  | ||  ||        |||| ||    |   || \n")
    cat("*********************")
    cat("*******************************************************\n")
    cat("Copyright (C) 2012-2017 Peter Carbonetto.\n")
    cat("See http://www.gnu.org/licenses/gpl.html for the full license.\n")
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
  # variances of the regression coefficients (s), and the posterior
  # mean estimates of the coefficients for the covariates (mu.cov),
  # which includes the intercept.
  logw   <- rep(0,ns)
  s      <- matrix(0,p,ns)
  mu.cov <- matrix(0,ncol(Z),ns)
                   
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
    out      <- outerloop(X,Z,y,family,SZy,SZX,c(sigma),c(sa),c(logodds),
                          c(alpha),c(mu),c(eta),tol,maxiter,verbose,NULL,
                          update.sigma,update.sa,optimize.eta,n0,sa0)
    logw     <- out$logw
    sigma    <- out$sigma
    sa       <- out$sa
    mu.cov[] <- out$mu.cov
    alpha[]  <- out$alpha
    mu[]     <- out$mu
    s[]      <- out$s
    eta[]    <- out$eta
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
        out <- outerloop(X,Z,y,family,SZy,SZX,sigma[i],sa[i],logodds[,i],
                         alpha[,i],mu[,i],eta[,i],tol,maxiter,verbose,i,
                         update.sigma,update.sa,optimize.eta,n0,sa0)
        logw[i]    <- out$logw
        sigma[i]   <- out$sigma
        sa[i]      <- out$sa
        mu.cov[,i] <- out$mu.cov
        alpha[,i]  <- out$alpha
        mu[,i]     <- out$mu
        s[,i]      <- out$s
        eta[,i]    <- out$eta
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
      if (optimize.eta)
        eta <- rep.col(eta[,i],ns)
      if (update.sigma)
        sigma <- rep(sigma[i],ns)
      if (update.sa)
        sa <- rep(sa[i],ns)
    }

    # Compute a variational approximation to the posterior distribution
    # for each candidate setting of the hyperparameters.
    if (verbose) {
      cat("Computing marginal likelihood for",ns,"combinations of",
          "hyperparameters.\n")
      cat("-iteration-   variational    max.   incl variance params\n")
      cat("outer inner   lower bound  change   vars   sigma      sa\n")
    }
    
    # For each setting of the hyperparameters, find a set of
    # parameters that locally minimize the K-L divergence between the
    # approximating distribution and the exact posterior.
    for (i in 1:ns) {
      out <- outerloop(X,Z,y,family,SZy,SZX,sigma[i],sa[i],logodds[,i],
                       alpha[,i],mu[,i],eta[,i],tol,maxiter,verbose,i,
                       update.sigma,update.sa,optimize.eta,n0,sa0)
      logw[i]    <- out$logw
      sigma[i]   <- out$sigma
      sa[i]      <- out$sa
      mu.cov[,i] <- out$mu.cov
      alpha[,i]  <- out$alpha
      mu[,i]     <- out$mu
      s[,i]      <- out$s
      eta[,i]    <- out$eta
    }
    if (verbose)
      cat("\n")
  }
    
  # (6) CREATE FINAL OUTPUT
  # -----------------------
  if (family == "gaussian") {
    fit <- list(family = family,n = n,n0 = n0,sa0 = sa0,mu.cov = mu.cov,
                update.sigma = update.sigma,update.sa = update.sa,
                prior.same = prior.same,optimize.eta = FALSE,logw = logw,
                sigma = sigma,sa = sa,logodds = logodds,alpha = alpha,
                mu = mu,s = s,eta = NULL)

    # Compute the proportion of variance in Y, after removing linear
    # effects of covariates, explained by the regression model.
    if (verbose)
      cat("Estimating proportion of variance in Y explained by model.\n");
    fit$model.pve <- varbvspve(X,fit,nr)

    # Compute the proportion of variance in Y, after removing linear
    # effects of covariates, explained by each variable.
    fit$pve           <- matrix(0,p,ns)
    rownames(fit$pve) <- colnames(X)
    sx                <- var1.cols(X)
    for (i in 1:ns) 
      fit$pve[,i] <- sx*(mu[,i]^2 + s[,i])/var1(y)
  } else if (family == "binomial")
    fit <- list(family = family,n = n,n0 = n0,mu.cov = mu.cov,sa0 = sa0,
                update.sigma = FALSE,update.sa = update.sa,
                optimize.eta = optimize.eta,prior.same = prior.same,
                logw = logw,sigma = NULL,sa = sa,logodds = logodds,
                alpha = alpha,mu = mu,s = s,eta = eta,model.pve = NA)
  
  # Add column names to some of the outputs.
  rownames(fit$alpha) <- colnames(X)
  rownames(fit$mu)    <- colnames(X)
  rownames(fit$s)     <- colnames(X)
  if (prior.same)
    fit$logodds <- c(fit$logodds)
  else
    rownames(fit$logodds) <- colnames(X)

  # Declare the return value as an instance of class 'varbvs'.
  class(fit) <- c("varbvs","list")
  return(fit)
}

# ----------------------------------------------------------------------
# This function implements one iteration of the "outer loop".
outerloop <- function (X, Z, y, family, SZy, SZX, sigma, sa, logodds,
                       alpha, mu, eta, tol, maxiter, verbose, outer.iter,
                       update.sigma, update.sa, optimize.eta, n0, sa0) {
  p <- ncol(X)
  if (length(logodds) == 1)
    logodds <- rep(logodds,p)

  # Note that we need to multiply the prior log-odds by log(10),
  # because varbvsnorm, varbvsbin and varbvsbinz define the prior
  # log-odds using the natural logarithm (base e).
  if (family == "gaussian") {

    # Optimize the variational lower bound for the Bayesian variable
    # selection model.
    out <- varbvsnorm(X,y,sigma,sa,log(10)*logodds,alpha,mu,tol,maxiter,
                      verbose,outer.iter,update.sigma,update.sa,n0,sa0)
    out$eta <- eta

    # Adjust the variational lower bound to account for integral over
    # the regression coefficients corresponding to the covariates.
    out$logw <-
      out$logw - determinant(crossprod(Z),logarithm = TRUE)$modulus/2
    
    # Compute the posterior mean estimate of the regression
    # coefficients for the covariates under the current variational
    # approximation.
    out$mu.cov <- c(with(out,SZy - SZX %*% (alpha*mu)))
  } else if (family == "binomial") {

    # Optimize the variational lower bound for the Bayesian variable
    # selection model.
    if (ncol(Z) == 1)
      out <- varbvsbin(X,y,sa,log(10)*logodds,alpha,mu,eta,tol,maxiter,
                       verbose,outer.iter,update.sa,optimize.eta,n0,sa0)
    else
      out <- varbvsbinz(X,Z,y,sa,log(10)*logodds,alpha,mu,eta,tol,maxiter,
                        verbose,outer.iter,update.sa,optimize.eta,n0,sa0)
    out$sigma <- sigma

    # Compute the posterior mean estimate of the regression
    # coefficients for the covariates under the current variational
    # approximation.
    Xr         <- with(out,c(X %*% (alpha*mu)))
    d          <- slope(out$eta)
    out$mu.cov <- c(solve(t(Z) %*% (Z*d),c((y - 0.5 - d*Xr) %*% Z)))
  }
  numiter  <- length(out$logw)
  out$logw <- out$logw[numiter]
  return(out)
}
