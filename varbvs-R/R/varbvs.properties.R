# Return the number of observations used to fit a model.
nobs.varbvs <- function (object, ...)
  object$n

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

# Return the names of the candidate variables.
labels.varbvs <- function (object, ...)
  rownames(object$alpha)

# Return the estimates of the regression coefficients at each
# hyperparameter setting, as well as the "averaged" estimates.
coef.varbvs <- function (object, ...) {
  ns <- length(object$w)
  if (ns == 1)
    out <- with(object,c(beta.cov,beta))
  else {
    out <- with(object,rbind(cbind(mu.cov,beta.cov),
                             cbind(mu,beta)))
    colnames(out)[ns+1] <- "Averaged"
  }
  return(out)
}
