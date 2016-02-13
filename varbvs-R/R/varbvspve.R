# Samples nr posterior estimates of the proportion of variance in Y
# explained by the Bayesian variable selection model fitted using a
# variational approximation. This function is only valid for the
# linear regression model (family = "gaussian") with an intercept.
varbvspve <- function (X, fit, nr = 1000) {

  # Get the number of variables (p) and the number of hyperparameter
  # settings (ns).
  p  <- ncol(X)
  ns <- length(fit$logw)

  # Check input X.
  if (!(is.matrix(X) & is.numeric(X) & sum(is.na(X)) == 0))
    stop("Input X must be a numeric matrix with no missing values.")
  if (nrow(fit$alpha) != p)
    stop("Inputs X and fit are not compatible.")
  
  # Initialize storage for posterior estimates of the proportion of
  # variance explained.
  pve <- rep(0,nr)

  # Compute the normalized (approximate) probabilities.
  w <- normalizelogweights(fit$logw)

  # For each sample, compute the proportion of variance explained.
  for (i in 1:nr) {

    # Draw a hyperparameter setting from the posterior distribution.
    j <- sample(ns,1,prob = w)
    
    # Sample the region coefficients.
    b <- with(fit,mu[,j] + sqrt(s[,j]) * rnorm(p))
    b <- b * (runif(p) < fit$alpha[,j])

    # Compute the proportion of variance explained.
    sz     <- c(var1(X %*% b))
    pve[i] <- sz/(sz + fit$sigma[j])
  }

  return(pve)
}
