# Select a subset of the candidate hyperparameter settings, and return
# a new varbvs object with these hyperparameter settings only.
subset.varbvs <- function (x, subset, ...) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(x,"varbvs"))
    stop("Input argument object must be an instance of class \"varbvs\".")

  # Get the unevaluated subset expression.
  subset <- substitute(subset)

  # Get the hyperparameter settings satisfying the 'subset' condition.
  i <- which(with(x,eval(subset)))

  # Output the new varbvs object with these hyperparameter settings
  # only.
  out         <- x
  out$sigma   <- out$sigma[i]
  out$sa      <- out$sa[i]
  out$logodds <- out$logodds[i]
  out$mu.cov  <- out$mu.cov[i]
  out$logw    <- out$logw[i]
  out$alpha   <- out$alpha[,i]
  out$mu      <- out$mu[,i]
  out$s       <- out$s[,i]
  out$pve     <- out$pve[,i]
  return(out)
}
