# Compute Bayes factors measuring improvement in "fit" when each
# candidate variable from "vars" is included in the model instead of
# variable i.
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
  mu           <- matrix(0,p,ns)
  s            <- matrix(0,p,ns)
  rownames(BF) <- rownames(fit$alpha)
  rownames(mu) <- rownames(fit$alpha)
  rownames(s)  <- rownames(fit$alpha)
  
  # COMPUTE PROXY PROBABILITIES
  # ---------------------------
  d <- diagsq(X)

  # Repeat for each hyperparameter setting.
  for (k in 1:ns) {

    # Get the hyperparameter values.
    sigma <- fit$sigma[k]
    sa    <- fit$sa[k]
    
    # Get the posterior means of the regression coefficients, and the
    # posterior inclusion probabilities.
    alpha0 <- fit$alpha[,k]
    mu0    <- fit$mu[,k]

    # Remove from the outcome y the linear effects of all variables
    # except for variable i.
    alpha0[i] <- 0
    y0        <- c(y - X %*% (alpha0*mu0))
    
    # Repeat for each candidate variable.
    for (j in vars) {

      # Add back in the linear effect of variable j (except in the
      # special case when i = j, because the effect was not removed in
      # the first place).
      if (i != j)
        y0 <- y0 + alpha0[j]*mu0[j] * X[,j]

      # Compute the Bayes factor.
      s[j,k]  <- sa*sigma/(sa*d[j] + 1)
      mu[j,k] <- s[j,k]*dot(y0,X[,j])/sigma
      BF[j,k] <- sqrt(s[j,k]/(sa*sigma)) * exp(mu[j,k]^2/(2*s[j,k]))

      # Remove the linear effect of variable j (except when i = j).
      if (i != j)
        y0 <- y0 - alpha0[j]*mu0[j] * X[,j]
    }
  }

  # Output the Bayes factors for the selected candidate proxy variables
  # only, as well as the estimated (posterior) means and variances of
  # the regression coefficients.
  return(list(BF = BF[vars,],mu = mu[vars,],s = s[vars,]))
}
