# Compute posterior statistics, ignoring correlations.
varbvsindep <- function (fit, X, Z, y) {

  # Get the number of samples (n), variables (p) and hyperparameter
  # settings (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(fit$logw)

  # Add intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)
  y <- c(as.double(y))
  
  # If necessary, convert the prior log-odds to a p x ns matrix.
  if (fit$prior.same)
    fit$logodds <- matrix(repmat(fit$logodds,p,1,byrow = TRUE))
                        
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior; see Chipman, George and
  # McCulloch, "The Practical Implementation of Bayesian Model
  # Selection," 2001.
  if (fit$family == "gaussian") {
    if (ncol(Z) == 1) {
      X <- X - rep.row(colMeans(X),n)
      y <- y - mean(y)
    } else {

      # The equivalent expressions in MATLAB are  
      #
      #   y = y - Z*((Z'*Z)\(Z'*y))
      #   X = X - Z*((Z'*Z)\(Z'*X))  
      #
      # This should give the same result as centering the columns of X
      # and subtracting the mean from y when we have only one
      # covariate, the intercept.
      y <- y - c(Z %*% solve(crossprod(Z),c(y %*% Z)))
      X <- X - Z %*% solve(crossprod(Z),t(Z) %*% X)
    }
  }

  # Initialize storage for the outputs.
  alpha <- matrix(0,p,ns)
  mu    <- zeros(0,p,ns)
  s     <- zeros(0,p,ns)

  # Calculate the mean (mu) and variance (s) of the coefficients given that
  # the coefficients are included in the model, and the posterior inclusion
  # probabilities (alpha), ignoring correlations between variables. Repeat
  # for each combination of the hyperparameters.
  for (i = 1:ns) {
    if (fit$family == "gaussian")
      # TO DO.
    else if (fit$family == "binomial")
      # TO DO.
  }  
}
