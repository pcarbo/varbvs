# For each variable j in "vars", this function returns a Bayes factor
# measuring the improvement in "fit" when variable j is included in
# the model instead of variable i; that is, a larger Bayes factor
# indicates a better model fit by swapping variables i and j. From an
# optimization perspective, this could be viewed as addressing the
# question, if you had to update the variational parameters for one
# variable so as to improve the "fit" of the variational
# approximation, after setting the posterior inclusion probability for
# variable i to zero, which variable would you choose?
#
# A few notes:
#
#   - The set of candidate variables "vars" could include variable i,
#     but does not have to.
#
#   - This only works for family = "gaussian".
#
#   - This function is currently not documented in the package.
varbvsproxybf <- function (X, Z, y, fit, i, vars) {

  # Get the number of samples (n), the number of variables (p), and
  # the number of hyperparameter settings (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(fit$w)

  # CHECK INPUTS
  # ------------
  # TO DO: Check all inputs.

  # Check input "fit".
  if (fit$family != "gaussian")
    stop(paste("Function varbvsproxybf is currently only implemented for",
               "linear regression (fit$family = \"gaussian\")"))
  
  # If the set of candidate variables is not provided, set it to all
  # the variables.
  if (missing(vars))
    vars <- 1:p
    
  # Add an intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)
  
  # PREPROCESSING STEPS
  # -------------------
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior.
  out <- remove.covariate.effects(X,Z,y)
  X   <- out$X
  y   <- out$y
  rm(out)
  
  # INITIALIZE STORAGE FOR OUTPUTS
  # ------------------------------
  # Create a p x ns matrix containing the Bayes factors.
  BF           <- matrix(0,p,ns)
  rownames(BF) <- rownames(fit$alpha)
  
  # COMPUTE PROXY PROBABILITIES
  # ---------------------------
  # Compute a couple useful quantities. 
  xy <- c(y %*% X)
  d  <- diagsq(X)

  # Repeat for each hyperparameter setting.
  for (k in 1:ns) {

    # Get the variational lower bound to the marginal log-likelihood.
    logw0 <- fit$logw[k]
      
    # Get the hyperparameter settings.
    sigma <- fit$sigma[k]
    sa    <- fit$sa[k]
    if (is.matrix(fit$logodds))
      logodds <- log(10)*fit$logodds[,k]
    else
      logodds <- rep(log(10)*fit$logodds[k],ncol(X))
    
    # Get the posterior mean and variances of the regression and the
    # posterior inclusion probabilities.
    alpha    <- fit$alpha[,k]
    mu       <- fit$mu[,k]
    s        <- fit$s[,k]

    # Set the probability for variable i to zero.
    alpha[i] <- 0
    
    # Repeat for each candidate variable.
    Xr <- c(X %*% (alpha * mu))
    for (j in vars) {

      # Save the current variational parameters for variable j.
      alpha0 <- alpha[j]
      mu0    <- mu[j]
        
      # Update the variational estimate of the posterior mean.
      r     <- alpha[j] * mu[j]
      mu[j] <- s[j]/sigma * (xy[j] + d[j]*r - sum(X[,j]*Xr))

      # Update the variational estimate of the posterior inclusion
      # probability.
      alpha[j] <- sigmoid(logodds[j] + (log(s[j]/(sa*sigma)) + mu[j]^2/s[j])/2)
      
      # Update Xr = X*r.
      Xr1 <- Xr + (alpha[j]*mu[j] - r) * X[,j]
      
      # Compute the variational estimate of the Bayes factor.
      #
      # TO DO: Speed up the calculation of the Bayes factor by
      # removing all the terms that automatically cancel out because
      # they are the unchanged after the update above.
      #
      logw1   <- int.linear(Xr1,d,y,sigma,alpha,mu,s) +
                 int.gamma(logodds,alpha) +
                 int.klbeta(alpha,mu,s,sigma*sa) -
                 determinant(crossprod(Z),logarithm = TRUE)$modulus/2
      BF[j,k] <- exp(logw1 - logw0)

      # Restore the variational parameters for variable j.
      alpha[j] <- alpha0
      mu[j]    <- mu0
    }
  }

  # Output the Bayes factors for the selected candidate variables only.
  return(BF[vars,])
}
