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
# Samples nr posterior estimates of the proportion of variance in Y
# explained by the Bayesian variable selection model fitted using a
# variational approximation. This function is only valid for the
# linear regression model (family = "gaussian") with an intercept.
varbvspve <- function (X, fit, nr = 1000) {

  # Get the number of variables (p) and the number of hyperparameter
  # settings (ns).
  p  <- ncol(X)
  ns <- length(fit$logw)

  # Check input X.
  if (!(is.matrix(X) & is.numeric(X) & sum(is.na(X)) == 0))
    stop("Input X must be a numeric matrix with no missing values.")
  if (nrow(fit$alpha) != p)
    stop("Inputs X and fit are not compatible.")

  # Check input "fit".
  if (!is(fit,"varbvs"))
    stop("Input argument \"fit\" must be an instance of class \"varbvs\".")
  if (fit$family != "gaussian")
    stop("varbvspve is only implemented for family = \"gaussian\"")
  
  # Initialize storage for posterior estimates of the proportion of
  # variance explained.
  pve <- rep(0,nr)

  # For each sample, compute the proportion of variance explained.
  for (i in 1:nr) {

    # Draw a hyperparameter setting from the posterior distribution.
    j <- sample(ns,1,prob = fit$w)
    
    # Sample the region coefficients.
    b <- with(fit,mu[,j] + sqrt(s[,j]) * rnorm(p))
    b <- b * (runif(p) < fit$alpha[,j])

    # Compute the proportion of variance explained.
    sz     <- c(var1(X %*% b))
    pve[i] <- sz/(sz + fit$sigma[j])
  }

  return(pve)
}
