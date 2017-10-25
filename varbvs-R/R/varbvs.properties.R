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
# Return the number of observations used to fit a model.
nobs.varbvs <- function (object, ...)
  nrow(object$fitted.values)

# ----------------------------------------------------------------------
# Return the names of the samples.
case.names.varbvs <- function (object, ...)
  names(object$fitted.values)

# ----------------------------------------------------------------------
# Return the names of the (included) variables.
variable.names.varbvs <- function (object, full = FALSE,
                                   include.threshold = 0.01, ...) {
  if (full)
    return(with(object,c(rownames(mu.cov),rownames(alpha))))
  else {
    i <- which(object$pip >= include.threshold)
    return(with(object,c(rownames(mu.cov),rownames(alpha)[i])))
  }
}

# ----------------------------------------------------------------------
# Return the names of the candidate variables.
labels.varbvs <- function (object, ...)
  rownames(object$alpha)

# ----------------------------------------------------------------------
# Return the estimates of the regression coefficients at each
# hyperparameter setting, as well as the "averaged" estimates.
coef.varbvs <- function (object, ...) {
  ns <- length(object$w)
  if (ns == 1)
    out <- with(object,c(beta.cov,beta))
  else {
    out <- with(object,rbind(cbind(mu.cov,beta.cov),
                             cbind(mu,beta)))
    colnames(out) <- c(paste0("theta_",1:ns),"Averaged")
  }
  return(out)
}

# ----------------------------------------------------------------------
# Extract the fitted values.
fitted.varbvs <- function (object, ...)
  object$fitted.values
