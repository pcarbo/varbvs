# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2018, Peter Carbonetto
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
# Generate a four-part summary of the fitted Bayesian variable
# selection model.
summary.varbvs <- function (object, cred.int = 0.95, nv, pip.cutoff, ...) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(object,"varbvs"))
    stop("Input argument object must be an instance of class \"varbvs\".")

  # Determine the number of variables for which to show detailed
  # summaries.
  if (missing(nv) & missing(pip.cutoff))
    nv <- 5
  else if (!missing(nv) & !missing(pip.cutoff))
    stop("Input arguments \"nv\" and \"pip.cutoff\" cannot both be specified")
  else if (!missing(pip.cutoff)) {
    if (pip.cutoff < 0 | pip.cutoff > 1)
      stop("Argument \"pip.cutoff\" should be a number between 0 and 1")
    nv <- sum(object$pip >= pip.cutoff)
  }
  
  # Get the normalized importance weights.
  w <- object$w
  
  # Get the number of variables (p) and number of candidate
  # hyperparameter settings (ns).
  p  <- nrow(object$alpha)
  ns <- length(w)

  # Input nv must be at least 1, and cannot be greater than the number
  # of variables.
  nv <- max(1,nv)
  nv <- min(p,nv) 
  
  # Generate the summary.
  out <-
    list(family       = object$family,
         cred.int     = cred.int,
         n            = nobs(object),
         p            = p,
         ns           = ns,
         ncov         = nrow(object$mu.cov),
         prior.same   = object$prior.same,
         update.sigma = object$update.sigma,
         update.sa    = object$update.sa,
         optimize.eta = object$optimize.eta,
         logw         = object$logw,
         w            = w,
         sigma        = list(x = NA,x0 = NA,a = NA,b = NA),
         sa           = list(x = NA,x0 = NA,a = NA,b = NA),
         logodds      = list(x = NA,x0 = NA,a = NA,b = NA),
         model.pve    = NULL)

  # Summarize the proportion of variance explained (PVE) by the fitted
  # model.
  if (!is.null(object$model.pve))
    out$model.pve =
      list(x0 = mean(object$model.pve),
           a  = quantile(object$model.pve,0.5 - cred.int/2,na.rm = TRUE),
           b  = quantile(object$model.pve,0.5 + cred.int/2,na.rm = TRUE))

  # Summarize the candidate hyperparameter settings, when provided.
  if (!object$update.sigma)
    out$sigma$x <- object$sigma
  if (!object$update.sa)
    out$sa$x <- object$sa
  if (object$prior.same)
    out$logodds$x <- object$logodds
  
  if (ns == 1) {

    # Summarize the hyperparameter settings when there is only one
    # candidate setting.
    out$sa$x0 <- object$sa
    if (object$family == "gaussian")
      out$sigma$x0 <- object$sigma
    if (object$prior.same)
      out$logodds$x0 <- object$logodds
  } else {
    
    # Summarize the residual variance parameter (sigma).
    if (object$family == "gaussian") {
      x <- object$sigma
      if (length(unique(x)) > 1) {
        x0        <- dot(w,x)
        out$sigma <- c(list(x = out$sigma$x,x0 = x0),cred(x,x0,w,cred.int))
      }
    }
 
    # Summarize the fitted prior variance parameter (sa).
    x <- object$sa
    if (length(unique(x)) > 1) {
      x0     <- dot(w,x)
      out$sa <- c(list(x = out$sa$x,x0 = x0),cred(x,x0,w,cred.int))
    }

    # Summarize the fitted prior log-odds of inclusion (logodds).
    if (object$prior.same) {
      x           <- object$logodds
      x0          <- dot(w,x)
      out$logodds <- c(list(x = out$logodds$x,x0 = x0),cred(x,x0,w,cred.int))
    }
  }

  # Summarize the number of variables selected at different PIP thresholds.
  out$num.included <-
      with(object,as.table(c(sum(pip>0.1),sum(pip>0.25),sum(pip>0.5),
                             sum(pip>0.75),sum(pip>0.9),sum(pip>0.95))))
  names(out$num.included) <- c(">0.10",">0.25",">0.50",">0.75",">0.90",">0.95")
  
  # Get more detailed statistics about the top nv variables by the
  # probability that they are included.
  vars <- order(object$pip,decreasing = TRUE)[1:nv]
  CIs  <- confint(object,vars,cred.int)
  CIs  <- do.call(rbind,lapply(CIs,function (x) x[ns + 1,]))
  out$top.vars <-
    data.frame(index = vars,variable = rownames(object$alpha)[vars],
               prob = object$pip[vars],PVE = NA,coef = object$beta[vars],
               cred = NA,stringsAsFactors = FALSE)
  for (i in 1:length(vars)) {
    if (!is.null(object$pve))
      out$top.vars[i,"PVE"] <- dot(w,object$pve[vars[i],])
    out$top.vars[i,"cred"] <- sprintf("[%+0.3f,%+0.3f]",CIs[i,1],CIs[i,2])
  }
  colnames(out$top.vars)[6] <- sprintf("Pr(coef.>%0.2f)",cred.int)
  rownames(out$top.vars)    <- NULL
  
  class(out) <- c("summary.varbvs","list")
  return(out)
}

# ----------------------------------------------------------------------
print.summary.varbvs <- function (x, digits = 3, ...) {

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
    if (!is.null(model.pve)) {
      with(model.pve,{
        cat("proportion of variance explained: ")
        cat(sprintf("%0.3f [%0.3f,%0.3f]\n",x0,a,b))
      })
    }
  
    # SUMMARIZE RESULTS ON HYPERPARAMETERS
    # ------------------------------------
    cat("Hyperparameters: ")
    if (ns == 1) {

      # Summarize the hyperparameter settings when there is only one
      # candidate setting.
      cat(sprintf("sigma=%0.3g sa=%0.3g ",sigma$x0,sa$x0))
      if (prior.same)
        cat(sprintf("logodds=%+0.2f",logodds$x0))
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
             cat(sprintf("logodds*%+8.2f %-19s (%+0.2f)--(%+0.2f)\n",x0,
                         sprintf("[%+0.2f,%+0.2f]",a,b),min(x),max(x))))
    }
    cat("*See help(varbvs) for details on how to convert between the\n")
    cat("prior log-odds and the prior inclusion probability.\n")
    
    # SUMMARIZE VARIABLE SELECTION RESULTS
    # ------------------------------------
    # Summarize the number of variables selected at different PIP thresholds.
    cat("Selected variables by probability cutoff:\n")
    print(num.included)
      
    # Give more detailed statistics about the top n variables by the
    # probability that they are included.
    if (x$family == "binomial")
      names(top.vars)[5] <- "coef*"
    cat(sprintf("Top %d variables by inclusion probability:\n",nrow(top.vars)))
    print(top.vars,digits = digits)
    if (x$family == "binomial")
      cat("*See help(varbvs) about interpreting coefficients in logistic",
          "regression.\n")
  })
  
  return(invisible(x))
}

# ----------------------------------------------------------------------
# Display summary of pcaviz object.
print.varbvs <- function (x, digits = 3, ...)
  print(summary(x),digits,...)
