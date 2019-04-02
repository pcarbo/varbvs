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
predict.varbvs <-
  function (object, X, Z = NULL, type = c("link","response","class"),
            averaged = TRUE, ...) {
  
  # Check that the first input is an instance of class "varbvs".
  if (!is(object,"varbvs"))
    stop("Input argument object must be an instance of class \"varbvs\".")

  # Process and check input argument "type".
  type <- match.arg(type)
  if (object$family == "gaussian" & type != "link")
    stop(paste("Prediction types \"response\" and \"class\" apply only to",
               "logistic regression (family = \"binomial\")"))
  
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
  if (ncol(Z) != nrow(object$mu.cov))
    stop("Inputs arguments object and Z are not compatible")

  # Get the normalized (approximate) probabilities.
  w <- object$w

  # Compute the estimates for each hyperparameter setting.
  out <- with(object,varbvs.linear.predictors(X,Z,mu.cov,alpha,mu))
  if (type == "response")
    out <- sigmoid(out)
  else if (type == "class")
    out <- round(sigmoid(out))
  
  # Average the estimates of Y over the hyperparameter settings, if
  # requested. For the logistic regression, the final "averaged"
  # prediction is obtained by collecting the "votes" from each
  # hyperparameter setting, weighting the votes by the marginal
  # probabilities, and outputing the estimate that wins by
  # majority. The averaged estimate is computed this way because the
  # estimates conditioned on each hyperparameter setting are not
  # necessarily calibrated in the same way.
  if (averaged) {
    if (type == "link")
      out <- c(out %*% w)
    else if (type == "response")
      out <- c(out %*% w)
    else
      out <- c(round(out %*% w))
  }
  names(out) <- rownames(X)
  return(out)
}

# ----------------------------------------------------------------------
# For each hyperparameter setting, and for each sample, compute a
# posterior mean estimate of Y. (For the logistic regression model, Y
# contains the posterior probability that the binary outcome is 1.)
varbvs.linear.predictors <- function (X, Z, mu.cov, alpha, mu) {
  ns <- ncol(alpha)
  Y  <- Z %*% mu.cov + X %*% (alpha*mu)
  return(Y)
}
