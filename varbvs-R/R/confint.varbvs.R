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
# Return "credible" or "confidence" intervals for all hyperparameter
# settings, as well as intervals averaging over all hyperparameter
# settings.
confint.varbvs <- function (object, parm, level = 0.95, ...) {

  # If input "parm" is not provided, select the top 5 variables by
  # posterior inclusion probability.
  if (missing(vars))
    vars <- order(object$pip,decreasing = TRUE)[1:5]
    
  # Get the number of hyperparameter settings (ns), the number of
  # selected variables (n), and the total number of variables (p).
  ns <- length(object$w)
  n  <- length(parm)
  p  <- nrow(object$alpha)
  
  # Check and process input "parm".
  variable.names <- rownames(object$alpha)
  if (is.numeric(parm))
    parm <- variable.names[parm]
  if (any(is.na(parm)) | !all(is.element(parm,variable.names)))
    stop(paste("Argument \"parm\" should contain only valid variable names",
               "or variable indices (columns of X)"))
  
  # Set up the data structure for storing the output.
  out        <- vector("list",n)
  names(out) <- variable.names
  
  # Compute the confidence intervals for for each requested variable.
  for (i in parm)
    out[[i]] <- get.confint.matrix(object,i,level)

  # No need to return a list if only one parameter (i.e., variable)
  # was requested. And in the special case when there is only 1
  # hyperparameter setting, return the confidence intervals in a
  # matrix.
  if (n == 1)
    out <- unlist(out,recursive = FALSE)
  else if (ns == 1)
    out <- do.call(rbind,out)
  return(out)
}

# Return an n x 2 matrix containing the x% confidence intervals (x =
# level) for variable i at each of the n hyperparameter settings, plus
# the confidence interval averaging over all hyperparameter settings.
get.confint.matrix <- function (fit, i, level) {

  # Get the number of hyperparameter settings.
  ns <- length(fit$w)
  
  if (ns == 1) {

    # In the special case when there is only one hyperparmaeter
    # setting, return the confidence interval in a 1 x 2 matrix.
    out <- credintnorm(level,fit$mu[i,],fit$s[i,])
    out <- matrix(out,1,2)
  } else {
  
    # Set up the data structure for storing the output.
    out           <- matrix(0,ns + 1,2)
    rownames(out) <- c(paste0("theta_",1:ns),"averaged")
        
    # Compute the credible interval for each hyperparameter setting.
    for (j in 1:ns)
      out[j,] <- credintnorm(level,fit$mu[i,j],fit$s[i,j])

    # Compute the credible interval averaging over all hyperparameter
    # settings.
    out[ns + 1,] <- credintmix(x,fit$w,fit$mu[i,],fit$s[i,])
  }

  # Add column labels indicating the lower and upper confidence limits.
  colnames(out) <- paste(round(100*c(0.5 - level/2,0.5 + level/2),
                               digits = 3),"%")
  return(out)
}

# Compute Monte Carlo estimates of credible intervals for coefficients
# in the fitted variable selection model. This function is used by
# summary.varbvs to generate credible intervals for coefficients of
# top-ranked variables.
varbvscoefcred <- function (fit, vars, cred.int = 0.95, nr = 1000) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(fit,"varbvs"))
    stop("Input \"fit\" must be an instance of class \"varbvs\"")
  
  # Get the number of hyperparameter settings.
  ns <- length(fit$logw)
  
  # Take care of optional inputs.
  if (missing(vars)) {
    p    <- nrow(fit$alpha)
    vars <- 1:p
  } else
    p <- length(vars)

  # Initialize storage for the result.
  a <- rep(0,p)
  b <- rep(0,p)
  
  # Repeat for each selected variable.
  for (i in 1:p) {
    j    <- vars[i]
    k    <- sample(ns,nr,prob = fit$w,replace = TRUE)
    x    <- fit$mu[j,k] + sqrt(fit$s[j,k]) * rnorm(nr)
    a[i] <- quantile(x,0.5 - cred.int/2)
    b[i] <- quantile(x,0.5 + cred.int/2)
  }
  
  return(list(a = a,b = b))
}
