# Generate a four-part summary of the fitted Bayesian variable
# selection model.
summary.varbvs <- function (fit, cred.int = 0.95, nv = 5, nr = 1000, ...) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(fit,"varbvs"))
    stop("Input fit must be an instance of class \"varbvs\".")
  
  # Get the number of variables (p) and number of candidate
  # hyperparameter settings (ns).
  p  <- nrow(fit$alpha)
  ns <- length(fit$logw)

  # Compute the normalized (approximate) probabilities.
  w <- c(normalizelogweights(fit$logw))

  # Compute the posterior inclusion probabilities (PIPs) and posterior
  # mean regression coefficients averaged over settings of the
  # hyperparameters.
  PIP  <- fit$alpha %*% w
  beta <- fit$mu    %*% w

  # Generate the summary.
  out <- list(family       = fit$family,
              cred.int     = cred.int,
              n            = fit$n,
              p            = p,
              ns           = ns,
              ncov         = fit$ncov,
              prior.same   = fit$prior.same,
              update.sigma = fit$update.sigma,
              update.sa    = fit$update.sa,
              optimize.eta = fit$optimize.eta,
              model.pve    = list(x0=mean(fit$model.pve),
                                  a=quantile(fit$model.pve,0.5 - cred.int/2),
                                  b=quantile(fit$model.pve,0.5 + cred.int/2)),
              sigma        = list(x = NA,x0 = NA,a = NA,b = NA),
              sa           = list(x = NA,x0 = NA,a = NA,b = NA),
              logodds      = list(x = NA,x0 = NA,a = NA,b = NA),
              logw         = fit$logw,
              w            = w)

  # Summarize the candidate hyperparameter settings, when provided.
  if (!fit$update.sigma)
    out$sigma$x <- fit$sigma
  if (!fit$update.sa)
    out$sa$x <- fit$sa
  if (fit$prior.same)
    out$logodds$x <- fit$logodds
  
  if (ns == 1) {

    # Summarize the hyperparameter settings when there is only one
    # candidate setting.
    out$sa$x0 <- fit$sa
    if (fit$family == "gaussian")
      out$sigma$x0 <- fit$sigma
    if (fit$prior.same)
      out$logodds$x0 <- fit$logodds
  } else {
    
    # Summarize the residual variance parameter (sigma).
    if (fit$family == "gaussian") {
      x <- fit$sigma
      if (length(unique(x)) > 1) {
        x0        <- dot(w,x)
        out$sigma <- c(list(x = out$sigma$x,x0 = x0),cred(x,x0,w,cred.int))
      }
    }
 
    # Summarize the fitted prior variance parameter (sa).
    x <- fit$sa
    if (length(unique(x)) > 1) {
      x0     <- dot(w,x)
      out$sa <- c(list(x = out$sa$x,x0 = x0),cred(x,x0,w,cred.int))
    }

    # Summarize the fitted prior log-odds of inclusion (logodds).
    if (fit$prior.same) {
      x           <- fit$logodds
      x0          <- dot(w,x)
      out$logodds <- c(list(x = out$logodds$x,x0 = x0),cred(x,x0,w,cred.int))
    }
  }

  # Summarize the number of variables selected at different PIP thresholds.
  out$n.incl <- as.table(c(sum(PIP > 0.1),sum(PIP > 0.25),sum(PIP > 0.5),
                           sum(PIP > 0.75),sum(PIP > 0.9),sum(PIP > 0.95)))
  names(out$n.incl) <- c(">0.10",">0.25",">0.50",">0.75",">0.90",">0.95")
  
  # Get more detailed statistics about the top nv variables by the
  # probability that they are included.
  vars         <- order(PIP,decreasing = TRUE)[1:nv]
  out$top.vars <-
    data.frame(index = vars,variable = rownames(fit$alpha)[vars],
               prob = PIP[vars],PVE = NA,coef = beta[vars],cred = NA)
  for (i in 1:length(vars)) {
    if (fit$family == "gaussian")
      out$top.vars[i,"PVE"] <- dot(w,fit$pve[vars[i],])
    out$top.vars[i,"cred"] <- with(varbvscoefcred(fit,vars[i],cred.int,nr),
                                   sprintf("[%+0.3f,%+0.3f]",a,b))
  }
  names(out$top.vars)[6] <- sprintf("Pr(coef.>%0.2f)",cred.int)
  
  class(out) <- c("summary.varbvs","list")
  return(out)
}

# ----------------------------------------------------------------------
print.summary.varbvs <- function (x, digits = 3) {

  # Check that the first input is an instance of class "summary.varbvs".
  if (!is(x,"summary.varbvs"))
    stop("Input must be an instance of class \"summary.varbvs\".")
  
  with(x,{
    
    # SUMMARIZE ANALYSIS SETUP
    # ------------------------
    cat("Summary of fitted Bayesian variable selection model:\n")
    cat(sprintf("family:     %-8s",family))
    cat(sprintf("   num. hyperparameter settings: %d\n",ns))
    cat(sprintf("samples:    %-6d",n))
    cat(sprintf("     iid variable selection prior: %s\n",tf2yn(prior.same)))
    cat(sprintf("variables:  %-6d",p))
    cat(sprintf("     fit prior var. of coefs (sa): %s\n",tf2yn(update.sa)))
    cat(sprintf("covariates: %-6d     ",ncov))
    if (family == "gaussian")
      cat(sprintf("fit residual var. (sigma):    %s\n",tf2yn(update.sigma)))
    else if (family == "binomial")
      cat(sprintf("fit approx. factors (eta):    %s\n",tf2yn(optimize.eta)))
    cat(sprintf("maximum log-likelihood lower bound: %0.4f\n",max(logw)))
    if (family == "gaussian") {
      with(model.pve,{
        cat("proportion of variance explained: ")
        cat(sprintf("%0.1f%% [%0.1f%%,%0.1f%%]\n",100*x0,100*a,100*b))
      })
    }
  
    # SUMMARIZE RESULTS ON HYPERPARAMETERS
    # ------------------------------------
    cat("Hyperparameters: ")
    if (ns == 1) {

      # Summarize the hyperparameter settings when there is only one
      # candidate setting.
      cat(sprintf("sigma=%0.3g sa=%0.3g",sigma,sa))
      if (prior.same)
        cat(sprintf("logodds=%+0.2f",logodds))
      cat("\n")
    } else {

      # Summarize the residual variance parameter (sigma).
      cat("\n")
      cat("        estimate ")
      cat(sprintf("Pr>%0.2f             candidate values\n",cred.int))
      if (family == "gaussian")
        with(sigma,cat(sprintf("sigma   %8.3g %-19s %0.3g--%0.3g\n",x0,
                               sprintf("[%0.3g,%0.3g]",a,b),min(x),max(x))))
 
      # Summarize the fitted prior variance parameter (sa).
      with(sa,cat(sprintf("sa      %8.3g %-19s %0.3g--%0.3g\n",x0,
                          sprintf("[%0.3g,%0.3g]",a,b),min(x),max(x))))
      
      # Summarize the fitted prior log-odds of inclusion (logodds).
      if (prior.same)
        with(logodds,
             cat(sprintf("logodds %+8.2f %-19s (%+0.2f)--(%+0.2f)\n",x0,
                         sprintf("[%+0.2f,%+0.2f]",a,b),min(x),max(x))))
    }

    # SUMMARIZE VARIABLE SELECTION RESULTS
    # ------------------------------------
    # Summarize the number of variables selected at different PIP thresholds.
    cat("Selected variables by probability cutoff:\n")
    print(n.incl)
      
    # Give more detailed statistics about the top n variables by the
    # probability that they are included.
    cat(sprintf("Top %d variables by inclusion probability:\n",nrow(top.vars)))
    print(top.vars,digits = digits)
  })
  
  return(invisible(fit))
}
