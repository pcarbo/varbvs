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
