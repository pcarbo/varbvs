intlinear <- function (Xr, d, y, sigma, alpha, mu, s) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # of the linear regression log-likelihood taken with respect to the
  # variational approximation. Inputs Xr and d must be equal to Xr =
  # X*(alpha*mu) and d = diag(X'*X). For a description of the
  # remaining inputs, see ‘varbvsoptimize’.
  n <- length(y)
  f <- -n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma) -
        dot(d,betavar(alpha,mu,s))/(2*sigma)
  return(f)
}

intklbeta <- function (alpha, mu, s, sa) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the negative
  # Kullback-Leibler divergence between the approximating distribution
  # and the prior of the coefficients. Input sa specifies the prior
  # variance of the coefficients. (It is not the same as the sa used
  # as input to ‘varbvsoptimize’.) See function ‘varbvsoptimize’ for
  # details on the inputs to this function.
  f <- (sum(alpha) + dot(alpha,log(s/sa)) - dot(alpha,s + mu^2)/sa)/2 -
        dot(alpha,log(alpha + eps)) -
        dot(1 - alpha,log(1 - alpha + eps))
  return(f)
}

intgamma <- function (logodds, alpha) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # on the prior inclusion probabilities taken with respect to the
  # variational approximation. The input arguments are: logodds, a
  # sqcalar, or a vector, specifying the prior log-odds; alpha,
  # mixture weights for the variational approximation. See function
  # ‘varbvsoptimize’ for details on the input arguments.

  # This is the same as 
  #
  #    sum(alpha*log(q) + (1-alpha)*log(1-q)).
  #  
  return(sum((alpha-1) * logodds + logsigmoid(logodds)))
}
