# TO DO: Explain what this function does, and how to use it.
#
# NOTE: The set of candidate variables could include variable i.
#
varbvsproxybf <- function (X, Z, y, fit, i, vars) {

  # CHECK INPUTS
  # ------------
  # TO DO.
  if (missing(vars))
    vars <- 1:ncol(X)
    
  # Get the number of samples (n), the number of variables (p), and
  # the number of hyperparameter settings (ns).
  n  <- length(y)
  p  <- ncol(X)
  ns <- length(fit$w)

  # Add intercept.
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
  out <- varbvs:::remove.covariate.effects(X,Z,y)
  X   <- out$X
  y   <- out$y
  rm(out)
  
  # INITIALIZE STORAGE FOR OUTPUTS
  # ------------------------------
  # TO DO: Explain what this matrix will contain.
  BF           <- matrix(0,p,ns)
  rownames(BF) <- rownames(fit$alpha)
  
  # COMPUTE PROXY PROBABILITIES
  # ---------------------------
  # Compute a couple useful quantities. 
  xy <- c(y %*% X)
  d  <- varbvs:::diagsq(X)

  # Repeat for each hyperparameter setting.
  for (k in 1:ns) {

    # Get the variational lower bound to the marginal log-likelihood.
    logw0 <- fit$logw[k]
      
    # Get the hyperparameter settings.
    sigma <- fit$sigma[k]
    sa    <- fit$sa[k]
    if (is.matrix(fit$logodds))
      logodds <- fit$logodds[,k]
    else
      logodds <- rep(fit$logodds[k],ncol(X))
    logodds <- log(10)*logodds
    
    # Get the posterior mean and variances of the regression and the
    # posterior inclusion probabilities, and set the probability for
    # variable i to zero.
    alpha    <- fit$alpha[,k]
    mu       <- fit$mu[,k]
    s        <- fit$s[,k]
    alpha[i] <- 0
    Xr       <- c(X %*% (alpha * mu))

    # Repeat for each candidate variable.
    for (j in vars) {

      # Save the current variational parameters for variable j.
      alpha0 <- alpha[j]
      mu0    <- mu[j]
        
      # Update the variational estimate of the posterior mean.
      r     <- alpha[j] * mu[j]
      mu[j] <- s[j]/sigma * (xy[j] + d[j]*r - sum(X[,j]*Xr))

      # Update the variational estimate of the posterior inclusion
      # probability.
      alpha[j] <-
        varbvs:::sigmoid(logodds[j] + (log(s[j]/(sa*sigma)) + mu[j]^2/s[j])/2)
      
      # Update Xr = X*r.
      Xrj <- Xr + (alpha[j]*mu[j] - r) * X[,j]
      
      # Compute the variational estimate of the Bayes factor.
      logw1   <- varbvs:::int.linear(Xrj,d,y,sigma,alpha,mu,s) +
                 varbvs:::int.gamma(logodds,alpha) +
                 varbvs:::int.klbeta(alpha,mu,s,sigma*sa) -
                 determinant(crossprod(Z),logarithm = TRUE)$modulus/2
      BF[j,k] <- exp(logw1 - logw0)

      # Restore the variational parameters for variable j.
      alpha[j] <- alpha0
      mu[j]    <- mu0
    }
  }

  # Output the Bayes factors.
  return(BF[vars,])
}
