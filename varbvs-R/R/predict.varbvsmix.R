# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2019, Peter Carbonetto
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
# Predict Y (outcome) given X (variables), Z (covariates) and model.
predict.varbvsmix <-
  function (object, X, Z = NULL, ...) {
  
  # Check that the first input is an instance of class "varbvsmix".
  if (!is(object,"varbvsmix"))
    stop("Input argument object must be an instance of class \"varbvsmix\".")

  # Get the number of samples (n), variables (p) and hyperparameter
  # settings (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(object$logw)
  
  # Check input X.
  if (!(is.matrix(X) & is.numeric(X) & sum(is.na(X)) == 0))
    stop("Input X must be a numeric matrix with no missing values.")
  storage.mode(X) <- "double"

  # Check input Z, and add an intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) | sum(is.na(Z)) > 0)
      stop("Input Z must be a numeric matrix with no missing values.")
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    storage.mode(Z) <- "double"
    Z <- cbind(1,Z)
  }
  if (ncol(Z) != length(object$mu.cov))
    stop("Inputs arguments object and Z are not compatible")

  # Compute the posterior mean estimates of the outcomes, Y.
  return(with(object,drop(Z %*% mu.cov + X %*% rowSums(alpha*mu))))
}
