# Compute fully-factorized variational approximation for Bayesian
# variable selection in linear (family = "gaussian") or logistic
# regression (family = "binomial"). See varbvs.Rd for details.
varbvs <- function (X, Z, y, family = "gaussian", sigma = NULL, sa = NULL,
                    update.sigma = NULL, update.sa = NULL, tol = 1e-4,
                    maxiter = 1e4, verbose = TRUE) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # (1) CHECK INPUTS
  # ----------------
  # If input Z is not NULL, it must have as many rows as X.
  if (!is.null(Z) & nrow(Z) != n)
      stop("Inputs X and Z do not match")

  # Add intercept.
  Z <- cbind(matrix(1,n,1),Z)

  # Input y must have n entries.
  if (length(y) != n)
    stop("Inputs X and y do not match")

  # Check choice of regression model.
  if (family != "gaussian" & family != "binomial")
    stop("family must be gaussian or binomial")

  # (2) PROCESS OPTIONS
  # -------------------
  # Get candidate settings for the variance of the residual (sigma),
  # if provided. Note that this option is not valid for a binary trait.
  if (is.null(sigma)) {
    sigma        <- var(y)
    update.sigma <- TRUE
  } else {
    if (is.null(update.sigma))
      update.sigma <- FALSE
    if (family == "binomial")
      stop("Input sigma is not allowed for family = binomial")
  }
  
  # Get candidate settings for the prior variance of the coefficients
  # (sa), if provided.
  if (is.null(sa)) {
    sa <- 1
    update.sa <- TRUE
  } else if (is.null(update.sa))
    update.sa <- FALSE
}
