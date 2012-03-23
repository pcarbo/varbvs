loginvgamma <- function (x, a, b) {
  # loginvgamma(x,a,b) returns the logarithm of the probability
  # density function of the inverse gamma distribution at elements of
  # x, with shape a and scale b.
  f <- a*log(b) - lgamma(a) - (a+1)*log(x+eps) - b/(x+eps)
  return(f)
}
  
logpve <- function (a, x) {
  # logpve(a,x) returns the logarithm of the density function for x,
  # given that r = a*x/(a*x + 1) is uniform on [0,1]. Input a must be
  # a positive scalar. This is useful for calculating the prior
  # probability of the variance parameter x, when we have a uniform
  # prior on r, the proportion of variance explained.
  r <- a*x/(a*x+1)          # Proportion of variance explained.
  f <- 2*log(r/x) - log(a)  # Log-density.
  return(f)
}
  
