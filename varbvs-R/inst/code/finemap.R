# TO DO: Explain what this function does, and how to use it.
#
# NOTE: The set of candidate variables could include variable i.
#
varbvsproxybf <- function (X, y, fit, i, vars) {

  # CHECK INPUTS
  # ------------
  # TO DO.
  if (missing(vars))
    vars <- 1:ncol(X)
    
  # Get the number of candidate variables (p) and the number of
  # hyperparameter settings (ns).
  p  <- length(vars)
  ns <- length(w)
    
  # Compute a couple useful quantities. 
  xy <- c(y %*% X)
  d  <- diagsq(X)

  # INITIALIZE STORAGE FOR OUTPUTS
  # ------------------------------
  # TO DO: Explain what this matrix will contain.
  BF           <- matrix(0,p,ns)
  rownames(BF) <- rownames(fit$alpha)[vars]
  
  # COMPUTE PROXY PROBABILITIES
  # ---------------------------
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
    
    # Get the posterior mean and variances of the regression and the
    # posterior inclusion probabilities, and set the probability for
    # variable i to zero.
    alpha    <- fit$alpha[,k]
    mu       <- fit$mu[,k]
    s        <- fit$s[,k]
    alpha[i] <- 0
    Xr       <- c(X %*% (alpha * mu))

    browser()

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
        varbvs::sigmoid(logodds[j] + (log(s[j]/(sa*sigma)) + mu[j]^2/s[j])/2)
      
      # Compute the variational estimate of the Bayes factor.
      logw1   <- int.linear(Xr,d,y,sigma,alpha,mu,s) +
                 int.gamma(logodds,alpha) +
                 int.klbeta(alpha,mu,s,sigma*sa)
      BF[j,k] <- exp(logw1 - logw0)

      # Restore the variational parameters for variable j.
      alpha[j] <- alpha0
      mu[j]    <- mu0
    }
  }
}
