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
# Compute posterior statistics, ignoring correlations.
varbvsindep <- function (fit, X, Z, y) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(fit,"varbvs"))
    stop("Input fit must be an instance of class \"varbvs\".")
  
  # Get the number of samples (n), variables (p) and hyperparameter
  # settings (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(fit$logw)

  # Check input X.
  if (!(is.matrix(X) & is.double(X) & sum(is.na(X)) == 0))
    stop("Input X must be a double-precision matrix with no missing values.")
  if (nrow(fit$alpha) != p)
    stop("Inputs X and fit are not compatible.")

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
  y <- c(as.double(y))
  
  # If necessary, convert the prior log-odds to a p x ns matrix.
  if (fit$prior.same)
    fit$logodds <- matrix(fit$logodds,p,ns,byrow = TRUE)
  
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior; see Chipman, George and
  # McCulloch, "The Practical Implementation of Bayesian Model
  # Selection," 2001.
  if (fit$family == "gaussian") {
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
      y <- y - c(Z %*% solve(crossprod(Z),c(y %*% Z)))
      X <- X - Z %*% solve(crossprod(Z),t(Z) %*% X)
    }
  }

  # Initialize storage for the outputs.
  alpha <- matrix(0,p,ns)
  mu    <- matrix(0,p,ns)
  s     <- matrix(0,p,ns)

  # Calculate the mean (mu) and variance (s) of the coefficients given that
  # the coefficients are included in the model, and the posterior inclusion
  # probabilities (alpha), ignoring correlations between variables. Repeat
  # for each combination of the hyperparameters.
  for (i in 1:ns) {
    if (fit$family == "gaussian")
      out <- with(fit,varbvsnormindep(X,y,sigma[i],sa[i],log(10)*logodds[,i]))
    else if (fit$family == "binomial")
      out <- with(fit,varbvsbinzindep(X,Z,y,eta[,i],sa[i],log(10)*logodds[,i]))
    alpha[,i] <- out$alpha
    mu[,i]    <- out$mu
    s[,i]     <- out$s
    rm(out)
  }

  return(list(alpha = alpha,mu = mu,s = s))
}

# ----------------------------------------------------------------------
# This function computes the mean (mu) and variance (s) of the
# coefficients given that they are included in the linear regression
# model, then it computes the posterior inclusion probabilities
# (alpha), ignoring correlations between variables. This function is
# used in 'varbvsindep', above.
varbvsnormindep <- function (X, y, sigma, sa, logodds) {
  s     <- sa*sigma/(sa*diagsq(X) + 1)
  mu    <- s*c(y %*% X)/sigma
  alpha <- sigmoid(logodds + (log(s/(sa*sigma)) + mu^2/s)/2)
  return(list(alpha = alpha,mu = mu,s = s))
}

# ----------------------------------------------------------------------
# This function computes the mean (mu) and variance (s) of the coefficients
# given that they are included in the logistic regression model, then it
# computes the posterior inclusion probabilities (alpha), ignoring
# correlations between variables. This function is used in varbvsindep.m.
varbvsbinzindep <- function (X, Z, y, eta, sa, logodds) {
  stats <- updatestats_varbvsbinz(X,Z,y,eta)
  s     <- sa/(sa*stats$xdx + 1)
  mu    <- s * stats$xy;
  alpha <- sigmoid(logodds + (log(s/sa) + mu^2/s)/2)
  return(list(alpha = alpha,mu = mu,s = s))
}
