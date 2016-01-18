# Generate a four-part summary of the fitted Bayesian variable
# selection model.
print.varbvs <- function (fit, cred.int = 0.95, n = 5, nr = 1000) {

  # Get the number of variables (p) and number of candidate
  # hyperparameter settings (ns).
  p  <- ncol(fit$alpha)
  ns <- length(fit$logw)

  # (1) COMPUTE POSTERIOR STATISTICS
  # --------------------------------
  # Compute the normalized (approximate) importance weights.
  w <- c(normalizelogweights(fit$logw))

  # Compute the posterior inclusion probabilities (PIPs) and posterior
  # mean regression coefficients averaged over settings of the
  # hyperparameters.
  PIP  <- fit$alpha %*% w
  beta <- fit$mu    %*% w

  # (2) SUMMARIZE ANALYSIS SETUP
  # ----------------------------
  cat("Summary of fitted Bayesian variable selection model:\n")
  cat(sprintf("family:     %-8s",fit$family))
  cat(sprintf("   num. hyperparameter settings: %d\n",length(fit$sa)))
  cat(sprintf("samples:    %-6d",fit$n))
  cat(sprintf("     iid variable selection prior: %s\n",tf2yn(fit$prior.same)))
  cat(sprintf("variables:  %-6d",p))
  cat(sprintf("     fit prior var. of coefs (sa): %s\n",tf2yn(fit$update.sa)))
  cat(sprintf("covariates: %-6d     ",fit$ncov))
  if (fit$family == "gaussian")
    cat(sprintf("fit residual var. (sigma):    %s\n",tf2yn(fit$update.sigma)))
  else if (fit$family == "binomial")
    cat(sprintf("fit approx. factors (eta):    %s\n",tf2yn(fit$optimize.eta)))
  cat(sprintf("maximum log-likelihood lower bound: %0.4f\n",max(fit$logw)))
  if (fit$family == "gaussian") {
    x  <- sort(fit$model.pve)
    x0 <- mean(x)
    a  <- x[floor((0.5 - cred.int/2)*length(x))]
    b  <- x[ceiling((0.5 + cred.int/2)*length(x))]
    cat("proportion of variance explained: ")
    cat(sprintf("%0.1f%% [%0.1f%%,%0.1f%%]\n",100*x0,100*a,100*b))
  }
  
  return(invisible(fit))
}
