# Predict Y (outcome) given X (variables), Z (covariates) and fitted model.
predict.varbvs <- function (fit, X, Z, ...) {
  
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
  if (ncol(Z) != nrow(fit$mu.cov))
    stop("Inputs Z and fit are not compatible")

  # Compute the normalized (approximate) probabilities.
  w <- c(normalizelogweights(fit$logw))
  
  # For each hyperparameter setting, and for each sample, compute the
  # posterior mean estimate of Y, and then average these estimates
  # over the hyperparameter settings. For the logistic regression, the
  # final "averaged" estimate is obtained by collecting the "votes"
  # from each hyperparameter setting, weighting the votes by the
  # marginal probabilities, and outputing the estimate that wins by
  # majority. The averaged estimate is computed this way because the
  # estimates conditioned on each hyperparameter setting are not
  # necessarily calibrated in the same way.
  Y <- with(fit,Z %*% mu.cov + X %*% (alpha*mu))
  if (fit$family == "gaussian")
    return(c(Y %*% w))
  else if (fit$family == "binomial")
    return(round(round(sigmoid(Y)) %*% w))
  else
    stop("Invalid setting for fit$family")
}
