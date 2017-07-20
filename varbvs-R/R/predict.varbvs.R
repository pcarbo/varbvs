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
# Predict Y (outcome) given X (variables), Z (covariates) and objectted model.
predict.varbvs <- function (object, X, Z = NULL, ...) {
  
  # Check that the first input is an instance of class "varbvs".
  if (!is(object,"varbvs"))
    stop("Input argument object must be an instance of class \"varbvs\".")

  # Get the number of samples (n), variables (p) and hyperparameter
  # settings (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(object$logw)
  
  # Check input X.
  if (!(is.matrix(X) & is.double(X) & sum(is.na(X)) == 0))
    stop("Input X must be a double-precision matrix with no missing values.")

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
  if (ncol(Z) != nrow(object$mu.cov))
    stop("Inputs arguments object and Z are not compatible")

  # Get the normalized (approximate) probabilities.
  w <- object$w
  
  # For each hyperparameter setting, and for each sample, compute the
  # posterior mean estimate of Y, and then average these estimates
  # over the hyperparameter settings. For the logistic regression, the
  # final "averaged" estimate is obtained by collecting the "votes"
  # from each hyperparameter setting, weighting the votes by the
  # marginal probabilities, and outputing the estimate that wins by
  # majority. The averaged estimate is computed this way because the
  # estimates conditioned on each hyperparameter setting are not
  # necessarily calibrated in the same way.
  Y <- with(object,Z %*% mu.cov + X %*% (alpha*mu))
  if (object$family == "gaussian")
    return(c(Y %*% w))
  else if (object$family == "binomial")
    return(round(c(round(sigmoid(Y)) %*% w)))
  else
    stop("Invalid setting for object$family")
}
