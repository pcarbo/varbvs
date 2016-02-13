# Compute Monte Carlo estimates of credible intervals for coefficients
# in the fitted variable selection model. This function is used by
# summary.varbvs to generate credible intervals for coefficients of
# top-ranked variables.
varbvscoefcred <- function (fit, vars, cred.int = 0.95, nr = 1000) {

  # Check that the first input is an instance of class "varbvs".
  if (!is(fit,"varbvs"))
    stop("Input fit must be an instance of class varbvs")
  
  # Get the number of hyperparameter settings.
  ns <- length(fit$logw)
  
  # Take care of optional inputs.
  if (missing(vars)) {
    p    <- nrow(fit$alpha)
    vars <- 1:p
  } else
    p <- length(vars)

  # Compute the normalized (approximate) probabilities.
  w <- normalizelogweights(fit$logw)

  # Initialize storage for the result.
  a <- rep(0,p)
  b <- rep(0,p)
  
  # Repeat for each selected variable.
  for (i in 1:p) {
    j    <- vars[i]
    k    <- sample(ns,nr,prob = w,replace = TRUE)
    x    <- fit$mu[j,k] + sqrt(fit$s[j,k]) * rnorm(nr)
    a[i] <- quantile(x,0.5 - cred.int/2)
    b[i] <- quantile(x,0.5 + cred.int/2)
  }
  
  return(list(a = a,b = b))
}
