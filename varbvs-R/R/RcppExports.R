sigmoid_rcpp <- function (x)
  .Call('varbvs_sigmoid_rcpp',x)

varbvsnormupdate_rcpp <- function(X, sigma, sa, logodds, xy, d,
                                  alpha, mu, Xr, i)
  invisible(.Call(C_varbvs_varbvsnormupdate_rcpp,X,sigma,sa,
                  logodds,xy,d,alpha,mu,Xr,i))
