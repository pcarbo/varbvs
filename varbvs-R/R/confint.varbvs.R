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
# Return "credible" or "confidence" intervals for all hyperparameter
# settings, as well as intervals averaging over all hyperparameter
# settings.
confint.varbvs <- function (object, parm, level = 0.95, ...) {

  # If input "parm" is not provided, select the top 5 variables by
  # posterior inclusion probability.
  if (missing(parm))
    parm <- order(object$pip,decreasing = TRUE)[1:5]
    
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
               "(column names of X) or variable indices (columns of X)"))
  
  # Set up the data structure for storing the output.
  out        <- vector("list",n)
  names(out) <- parm
  
  # Compute the confidence intervals for for each requested variable.
  for (i in parm)
    out[[i]] <- get.confint.matrix(object,i,level)
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
    # setting, return the confidence interval in a 2 x 2 matrix.
    out <- credintnorm(level,fit$mu[i,],fit$s[i,])
    out <- matrix(out,2,2,byrow = TRUE)
  } else {
  
    # Set up the data structure for storing the output.
    out           <- matrix(0,ns + 1,2)
        
    # Compute the credible interval for each hyperparameter setting.
    for (j in 1:ns)
      out[j,] <- credintnorm(level,fit$mu[i,j],fit$s[i,j])

    # Compute the credible interval averaging over all hyperparameter
    # settings.
    out[ns + 1,] <- credintmix(level,fit$w,fit$mu[i,],fit$s[i,])
  }

  # Add row labels, and column labels indicating the lower and upper
  # confidence limits.
  rownames(out) <- c(colnames(fit$alpha),"averaged")
  colnames(out) <- paste(round(100*c(0.5 - level/2,0.5 + level/2),
                               digits = 3),"%")
  return(out)
}
