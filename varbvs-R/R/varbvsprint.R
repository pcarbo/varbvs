# Generate a four-part summary of the fitted Bayesian variable
# selection model.
varbvsprint <- function (fit, cred.int = 0.95, n = 5, nr = 1000) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(fit,"varbvs"))
    stop("Input fit must be an instance of class varbvs")
  
  # Get the number of variables (p) and number of candidate
  # hyperparameter settings (ns).
  p  <- nrow(fit$alpha)
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
    x <- fit$model.pve
    cat("proportion of variance explained: ")
    cat(sprintf("%0.1f%% [%0.1f%%,%0.1f%%]\n",100*mean(x),
                100*quantile(x,0.5 - cred.int/2),
                100*quantile(x,0.5 + cred.int/2)))
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
    cat("        estimate ")
    cat(sprintf("Pr>%0.2f             candidate values\n",cred.int))
    if (fit$family == "gaussian") {
      x0 <- dot(w,fit$sigma)
      cat(sprintf("sigma   %8.3g ",x0))
      cat(sprintf("%-19s ",with(cred(fit$sigma,x0,w,cred.int),
                                sprintf("[%0.3g,%0.3g]",a,b))))
      if (fit$update.sigma)
        cat("NA\n")
      else
        with(fit,cat(sprintf("%0.3g--%0.3g\n",min(sigma),max(sigma))))
    }
 
    # Summarize the fitted prior variance parameter (sa).
    x0 <- dot(w,fit$sa)
    cat(sprintf("sa      %8.3g ",x0))
    cat(sprintf("%-19s ",with(cred(fit$sa,x0,w,cred.int),
                              sprintf("[%0.3g,%0.3g]",a,b))))
    if (fit$update.sa)
      cat("NA\n")
    else
      with(fit,cat(sprintf("%0.3g--%0.3g\n",min(sa),max(sa))))

    # Summarize the fitted prior log-odds of inclusion (logodds).
    if (fit$prior.same) {
      x  <- fit$logodds
      x0 <- dot(w,x)
      cat(sprintf("logodds %+8.2f %-19s (%+0.2f)--(%+0.2f)\n",x0,
                  with(cred(x,x0,w,cred.int),
                       sprintf("[%+0.2f,%+0.2f]",a,b)),min(x),max(x)))
    }
  }

  # (4) SUMMARIZE VARIABLE SELECTION RESULTS
  # ----------------------------------------
  # Summarize the number of variables selected at different PIP thresholds.
  cat("Selected variables:\n")
  cat("prob. >0.10 >0.25 >0.50 >0.75 >0.90 >0.95\n")
  cat(sprintf("count %5d %5d %5d %5d %5d %5d\n",sum(PIP > 0.1),sum(PIP > 0.25),
              sum(PIP > 0.5),sum(PIP > 0.75),sum(PIP > 0.9),sum(PIP > 0.95)))
      
  # Give more detailed statistics about the top n variables by the
  # probability that they are included.
  vars <- order(PIP,decreasing = TRUE)[1:n]
  cat(sprintf("Top %d variables by inclusion probability:\n",n))
  cat(" index variable   prob.")
  if (fit$family == "gaussian")
    cat(" -PVE-")
  cat(sprintf("   coef. Pr(coef.>%0.2f)\n",cred.int))
  for (i in vars) {
    cat(sprintf("%6d %-10s %0.3f",i,rownames(fit$alpha)[i],PIP[i]))
    if (fit$family == "gaussian")
      cat(sprintf(" %04.1f%%",100*dot(w,fit$pve[i,])))
    with(varbvscoefcred(fit,i,cred.int,nr),
         cat(sprintf(" %+7.3f [%+0.3f,%+0.3f]\n",beta[i],a,b)))
  }

  return(invisible(fit))
}
