# varbvsmix.m: Fit linear regression model with mixture-of-normals
# prior using variational approximation. See varbvsmix.Rd for details.
varbvsmix <- function (eX, Z, y, sa, sigma, q, alpha, mu, update.sigma,
                       update.sa, update.q, q.penalty, tol = 1e-4,
                       maxiter = 1e4, verbose = TRUE) {
    
  # Get the number of samples (n), variables (p) and mixture
  # components (K).
  n <- nrow(X)
  p <- ncol(X)
  K <- length(sa)

  # (1) CHECK INPUTS
  # ----------------
  # Check input X.
  if (!(is.matrix(X) & is.double(X) & sum(is.na(X)) == 0))
    stop("Input X must be a double-precision matrix with no missing values.")

  # Add column names to X if they are not provided.
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X",1:p)

  # Check input Z.
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) | sum(is.na(Z)) > 0)
      stop("Input Z must be a numeric matrix with no missing values.")
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    storage.mode(Z) <- "double"
  }

  # Add intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)

  # Check input y.
  if (!is.numeric(y) | sum(is.na(y)) > 0)
    stop("Input y must be a numeric vector with no missing values.")
  if (length(y) != n)
    stop("Inputs X and y do not match")
  y <- c(as.double(y))
  
  # Check input sa. The variance of the first mixture component should
  # be exactly zero, corresponding to the "spike".
  sa <- c(as.double(sa))
  if (sa[1] != 0)
    stop("Variance of first mixture component should be 0")
}
