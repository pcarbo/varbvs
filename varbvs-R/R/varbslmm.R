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
  # effects, if provided.
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
  
  # Determine whether to update the residual variance parameter.
  if (missing(update.sigma))
    update.sigma <- update.sigma.default

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
}
