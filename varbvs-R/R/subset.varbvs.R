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
# Select a subset of the candidate hyperparameter settings, and return
# a new varbvs object with these hyperparameter settings only.
subset.varbvs <- function (x, subset, ...) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(x,"varbvs"))
    stop("Input argument object must be an instance of class \"varbvs\".")

  # Get the unevaluated subset expression.
  e <- substitute(subset)

  # Get the hyperparameter settings satisfying the 'subset' condition.
  i <- which(eval(e,x,parent.frame()))
  if (length(i) == 0)
    stop("No hyperparameter settings are selected.")

  # Output the new varbvs object with these hyperparameter settings
  # only.
  out         <- x
  out$sa      <- out$sa[i]
  out$logodds <- out$logodds[i]
  out$logw    <- out$logw[i]
  out$w       <- normalizelogweights(out$logw)
  out$mu.cov  <- as.matrix(out$mu.cov[,i])
  out$alpha   <- as.matrix(out$alpha[,i])
  out$mu      <- as.matrix(out$mu[,i])
  out$s       <- as.matrix(out$s[,i])
  out$fitted.values <- as.matrix(out$fitted.values)
  out$residuals     <- as.matrix(out$residuals)
  if (!is.null(out$pve))
    out$pve <- as.matrix(out$pve[,i])
  if (out$family == "gaussian")
    out$sigma <- out$sigma[i]
  else if (out$family == "binomial")
    out$eta <- out$eta[,i]
  return(out)
}
