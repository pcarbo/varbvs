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
# Return the number of observations used to fit a model.
nobs.varbvs <- function (object, ...)
  nrow(object$fitted.values)

# ----------------------------------------------------------------------
# Return the names of the samples.
case.names.varbvs <- function (object, ...)
  rownames(object$fitted.values)

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
# hyperparameter setting, as well as the "averaged" estimates in the
# last column of this matrix.
coef.varbvs <- function (object, ...) {
  ns <- length(object$w)
  if (ns == 1)
    out <- with(object,rbind(mu.cov,alpha*mu))
  else {
    out <- with(object,rbind(cbind(mu.cov,beta.cov),
                             cbind(alpha*mu,beta)))
    colnames(out)[ns + 1] <- "averaged"
  }
  return(out)
}

# ----------------------------------------------------------------------
# Return the fitted values stored in an n x ns matrix, where n is the
# number of samples and ns is the number of hyperparameter settings.
fitted.varbvs <- function (object, ...)
  object$fitted.values

# ----------------------------------------------------------------------
# Return the residuals stored in an n x ns matrix, where n is the
# number of samples and ns is the number of hyperparameter
# settings. For a logistic regression model, there are two types of
# residuals ("deviance" and "response").
resid.varbvs <- function (object, type = c("deviance","response"), ...) {
  if (object$family == "gaussian")
    out <- object$residuals
  else {
    type <- match.arg(type)
    if (type == "deviance")
      out <- object$residuals$deviance
    else if (type == "response")
      out <- object$residuals$response
    else
      stop("Argument \"type\" should be \"deviance\" or \"response\"")
  }
  return(out)
}
residuals.varbvs <- resid.varbvs

# ----------------------------------------------------------------------
# Return the deviance for each hyperparameter setting.
deviance.varbvs <- function (object, ...)
  colSums(resid(object,type = "deviance")^2)
