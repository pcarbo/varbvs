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
  out$mu.cov  <- as.matrix(out$mu.cov[,i])
  out$alpha   <- as.matrix(out$alpha[,i])
  out$mu      <- as.matrix(out$mu[,i])
  out$s       <- as.matrix(out$s[,i])
  if (out$family == "gaussian") {
    out$sigma <- out$sigma[i]
    out$pve   <- as.matrix(out$pve[,i])
  } else if (out$family == "binomial")
    out$eta <- out$eva[,i]
  return(out)
}
