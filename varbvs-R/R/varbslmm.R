# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2017, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# TO DO: Add brief description of function here.
#
# NOTES:
#
#   - "a" in sa is for "additive effects".
#
#   - "b" in sa is for "background effects".
# 
#   - Need to treat special case when sa = 0.
#
varbslmm <- function (X, Z, y, sigma, sa, sb, logodds, alpha, mu,
                      update.sigma, update.sa, initialize.params,
                      sa0 = 1, n0 = 10, tol = 1e-4, maxiter = 1e4,
                      verbose = TRUE) {

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

  # (2) PROCESS OPTIONS
  # -------------------
  if (!is.finite(maxiter))
    stop("Input maxiter must be a finite number")
  
  # Get candidate settings for the variance of the residual (sigma),
  # if provided.
  if (missing(sigma)) {
    sigma <- var(y)
    update.sigma.default <- TRUE
  } else {
    sigma <- c(sigma)
    update.sigma.default <- FALSE
  }

  # Get candidate settings for the prior variance of the "background"
  # effects, if provided.
  if (missing(sb))
    sb <- 0.01
  sb <- c(sb)

  # Get candidate settings for the prior variance of the "large"
  # additive effects, if provided.
  if (missing(sa)) {
    sa <- 1
    update.sa.default <- TRUE
  } else {
    sa <- c(sa)
    update.sa.default <- FALSE
  }
  
  # Get candidate settings for the prior log-odds of inclusion. A default
  # setting is only available if the number of other hyperparameter
  # settings is 1, in which case we select 20 candidate settings for
  # the prior log-odds, evenly spaced between log10(1/p) and -1.
  if (missing(logodds)) {
    if (length(sigma) == 1 & length(sa) == 1 & length(sb) == 1)
      logodds <- seq(-log10(p),-1,length.out = 20)
    else
      stop(paste("logodds can only be missing when length(sigma) =",
                 "length(sa) = length(sb) = 1")
  }
  logodds <- c(logodds)
  
  # Here is where I ensure that the numbers of candidate hyperparameter
  # settings are consistent.
  ns <- max(length(sigma),length(logodds),length(sa),length(sb))
  if (length(sigma) == 1)
    sigma <- rep(sigma,ns)
  if (length(logodds) == 1)
    logodds <- rep(logodds,ns)
  if (length(sa) == 1)
    sa <- rep(sa,ns)
  if (length(sb) == 1)
    sb <- rep(sb,ns)
  if (length(sigma) != ns | length(logodds) != ns |
      length(sa) != ns | length(sb) != ns)
    stop("Arguments sigma, logodds, sa and sb are inconsistent")

  # Determine whether to update the residual variance parameter.
  if (missing(update.sigma))
    update.sigma <- update.sigma.default

  # Determine whether to update the prior variance of the "large"
  # additive effects.
  if (missing(update.sa))
    update.sa <- update.sa.default

  # Set initial estimates of variational parameter alpha.
  initialize.params.default <- TRUE
  if (missing(alpha)) {
    alpha <- rand(p,ns)
    alpha <- alpha / rep.row(colSums(alpha),p)
  } else
    initialize.params.default <- FALSE
  if (!is.matrix(alpha))
    alpha <- matrix(alpha)
  if (nrow(alpha) != p)
    stop("Input alpha must have as many rows as X has columns")
  if (ncol(alpha) == 1)
    alpha <- rep.col(alpha,ns)

  # Set initial estimates of variational parameter mu.
  if (missing(mu))
    mu <- randn(p,ns)
  else
    initialize.params.default <- FALSE    
  if (!is.matrix(mu))
    mu <- matrix(mu)
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
    
  # (3) PREPROCESSING STEPS
  # -----------------------
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior.
  out <- remove.covariate.effects(X,Z,y)
  X   <- out$X
  y   <- out$y
  SZy <- out$SZy
  SZX <- out$SZX
  rm(out)
  
  # Compute the kinship matrix.
  fprintf('Computing kinship matrix.\n');
  K = tcrossprod(X)/p

  # Provide a brief summary of the analysis.
  if (verbose) {
    cat("Fitting variational approximation for Bayesian sparse linear",
        "mixed model.\n")
    cat(sprintf("samples:      %-6d ",n))
    cat(sprintf("mixture component sd's:         %0.2g..%0.2g\n",
                min(sqrt(sa[2:K])),max(sqrt(sa[2:K]))))
    cat(sprintf("variables:    %-6d ",p))
    cat(sprintf("mixture component drop thresh.: %0.1e\n",drop.threshold))
    cat(sprintf("covariates:   %-6d ",max(0,ncol(Z) - 1)))
    cat(sprintf("fit mixture weights:            %s\n",tf2yn(update.w)))
    cat(sprintf("mixture size: %-6d ",K))
    cat(sprintf("fit residual var. (sigma):      %s\n",tf2yn(update.sigma)))
    cat("intercept:    yes    ")
    cat(sprintf("convergence tolerance      %0.1e\n",tol))
  }

  # (4) INITIALIZE STORAGE FOR THE OUTPUTS
  # --------------------------------------
  # TO DO.
    
  # (5) FIT BSLMM MODEL TO DATA
  # ---------------------------
  if (ns == 1) {

    # Find a set of parameters that locally minimize the K-L
    # divergence between the approximating distribution and the exact
    # posterior.
    if (verbose) {
      # TO DO.
    }
    # TO DO.
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
        # TO DO.
      }

      # Repeat for each setting of the hyperparameters.
      for (i in 1:ns) {
        # TO DO.
      }
      if (verbose)
        cat("\n")

      # Choose an initialization common to all the runs of the
      # coordinate ascent algorithm. This is chosen from the
      # hyperparameters with the highest variational estimate of the
      # marginal likelihood.
      # TO DO.
    }

    # Compute a variational approximation to the posterior distribution
    # for each candidate setting of the hyperparameters.
    if (verbose) {
      cat("Computing marginal likelihood for",ns,"combinations of",
          "hyperparameters.\n")
      # TO DO.
    }
    
    # For each setting of the hyperparameters, find a set of
    # parameters that locally minimize the K-L divergence between the
    # approximating distribution and the exact posterior.
    for (i in 1:ns) {
      # TO DO.
    }
    if (verbose)
      cat("\n")
  }
    
  # (6) CREATE FINAL OUTPUT
  # -----------------------
  # TO DO.
    
  browser()
}

# ----------------------------------------------------------------------
# This function implements one iteration of the "outer loop" for the
# varbslmm function.
varbvslmm.outerloop <-
  function (X, K, y, alpha, mu, sigma, logodds, sa, update.sigma, update.sa,
            tol, maxiter, verbose, outer.iter, n0, sa0) {

}
