# Implements the fully-factorized variational approximation for
# Bayesian variable selection in linear regression. It finds the
# "best" fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors
# on the coefficients. By "best", we mean the approximating
# distribution that locally minimizes the K-L divergence between the
# approximating distribution and the exact posterior.
#
# Input X is an n x p matrix of observations about the variables (or
# features), where n is the number of samples, and p is the number of
# variables. Input y is the vector of observations about the outcome;
# it is a vector of length n.
#
# Inputs sigma, sa and logodds are additional model parameters; sigma
# and sa are scalars. Input sigma specifies the variance of the
# residual, and sa specifies the prior variance of the
# coefficients. logodds is the prior log-odds of inclusion for each
# variable.
#
# Output logw is the variational estimate of the marginal
# log-likelihood given the hyperparameters. Outputs alpha, mu and s
# are the parameters of the variational approximation and,
# equivalently, variational estimates of posterior quantites: under
# the variational approximation, the ith regression coefficient is
# normal with probability alpha(i); mu(i) and s(i) are the mean and
# variance of the coefficient given that it is included in the model.
#
# When update.sa = TRUE, there is the additional option of computing
# the maximum a posteriori (MAP) estimate of the prior variance
# parameter (sa), in which sa is drawn from a scaled inverse
# chi-square distribution with scale sa0 and degrees of freedom n0.
varbvsnorm <-
  function varbvsnorm (X, y, sigma, sa, logodds, alpha, mu, tol, maxiter, 
                       verbose, outer.iter, update.sigma, update.sa, n0, sa0) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)
  
  # (1) INITIAL STEPS
  # -----------------
  # Compute a few useful quantities. 
  xy <- c(y %*% X)
  d  <- diagsq(X)
  Xr <- c(X %*% (alpha*mu))
  
  # Calculate the variance of the coefficients.
  s <- sa*sigma/(sa*d + 1)

  # (2) MAIN LOOP
  # -------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  logw <- (-Inf)
  for (iter in 1:maxiter) {

    # Save the current variational parameters.
    alpha0 <- alpha
    mu0    <- mu

    # (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logw0 <- int_linear(Xr,d,y,sigma,alpha,mu,s) ...
             + int_gamma(logodds,alpha) ...
             + int_klbeta(alpha,mu,s,sigma*sa);
  }
}
  


























