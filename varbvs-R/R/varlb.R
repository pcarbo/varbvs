intlinear <- function (Xr, d, y, sigma, alpha, mu, s) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # of the linear regression log-likelihood taken with respect to the
  # variational approximation. Inputs Xr and d must be equal to Xr =
  # X*(alpha*mu) and d = diag(X'*X). For a description of the
  # remaining inputs, see ‘varbvsoptimize’.
  n <- length(y)
  f <- (-n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma)
        - dot(d,betavar(alpha,mu,s))/(2*sigma))
  return(f)
}

intlinearmix <- function (Xr, d, y, sigma, alpha, mu1, mu2, s1, s2) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood for the Bayesian variable selection
  # mixture model. This integral is the expectation of the linear
  # regression log-likelihood taken with respect to the variational
  # approximation to the mixture model.
  n <- length(y)
  f <- (-n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma) 
        - dot(d,betavarmix(alpha,mu1,mu2,s1,s2))/(2*sigma))
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
  f <- ((sum(alpha) + dot(alpha,log(s/sa)) - dot(alpha,s + mu^2)/sa)/2
        - dot(alpha,log(alpha + eps))
        - dot(1 - alpha,log(1 - alpha + eps)))
  return(f)
}

intklbetamix <- function (alpha, mu1, mu2, s1, s2, sa1, sa2) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood for the mixture model. This
  # integral is the negative Kullback-Leibler divergence between the
  # approximating distribution and the prior of the coefficients.
  a <- alpha
  f <- ((sum(a) + dot(a,log(s1/sa1)) - dot(a,s1 + mu1^2)/sa1)/2 
        + (sum(1-a) + dot(1-a,log(s2/sa2)) - dot(1-a,s2 + mu2^2)/sa2)/2 
        - dot(a,log(a + eps)) - dot(1-a,log(1-a + eps)))
  
  return(f)
}
