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

  # (3) SUMMARIZE RESULTS ON HYPERPARAMETERS
  # ----------------------------------------
  # Summarize the fitted residual variance parameter (sigma).
  cat("Hyperparameters: ")
  if (ns == 1) {

    # Summarize the hyperparameter settings when there is only one
    # candidate setting.
    if (fit$family == "gaussian")
      cat(sprintf("sigma=%0.3g ",fit$sigma))
    cat(sprintf("sa=%0.3g ",fit$sa))
    if (fit$prior.same)
      cat(sprintf("logodds=%+0.2f",fit$logodds))
    cat("\n")
  } else {
    cat("\n")
    cat(sprintf("        estimate Pr>%0.2f             candidate values\n",
                cred.int))
    if (fit$family == "gaussian") {
      x0  <- dot(w,fit$sigma)
      out <- cred(fit$sigma,x0,w,cred.int)
      cat(sprintf("sigma   %8.3g %-19s ",x0,sprintf("[%0.3g,%0.3g]",
                                                    out$a,out$b)))
      if (fit$update.sigma)
        cat("NA\n")
      else
        cat(sprintf("%0.3g--%0.3g\n",min(fit$sigma),max(fit$sigma)))
    }
 
    # Summarize the fitted prior variance parameter (sa).
    x0 <- dot(w,fit$sa)
    # TO DO: FIX THIS.
    # [a b] = cred(fit$sa,x0,w,c)
    a <- 0
    b <- 0
    cat(sprintf("sa      %8.3g %-19s ",x0,sprintf("[%0.3g,%0.3g]",a,b)))
    if (fit$update.sa)
      cat("NA\n")
    else
      cat(sprintf("%0.3g--%0.3g\n",min(fit$sa),max(fit$sa)))

    # Summarize the fitted prior log-odds of inclusion (logodds).
    if (fit$prior.same) {
      x  <- fit$logodds
      x0 <- dot(w,x)
      # TO DO: FIX THIS.
      # [a b] = cred(x,x0,w,c)
      a <- 0
      b <- 0
      cat(sprintf("logodds %+8.2f %-19s (%+0.2f)--(%+0.2f)\n",x0,
                  sprintf("[%+0.2f,%+0.2f]",a,b),min(x),max(x)))
    }
  }
  
  return(invisible(fit))
}
