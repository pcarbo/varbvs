# Generate a four-part summary of the fitted Bayesian variable
# selection model.
print.varbvs <- function (fit, c = 0.95, n = 5, nr = 1000) {

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
 
  return(NULL)
}
