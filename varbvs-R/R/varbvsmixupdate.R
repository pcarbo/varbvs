# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2017, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# [alpha,mu,Xr] = varbvsmixupdate(X,sigma,sa,q,xy,d,alpha0,mu0,Xr0,i)
# runs a single iteration of the coordinate ascent updates maximizing the
# variational lower bound for the linear regression model with a
# mixture-of-normals prior.
varbvsmixupdate <- function (X, sigma, sa, w, xy, d, alpha0, mu0, Xr0, i) {

  # Get the number of samples (n), the number of variables (p), and the
  # number of mixture components including the "spike" (K).
  n <- nrow(X)
  p <- ncol(X)
  K <- length(w)

  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input X should be a double-precision matrix")

  # Check input sigma.
  if (length(sigma) != 1)
    stop("Input sigma should be a scalar")

  # Check input sa.
  if (length(sa) != K)
    stop("Input sa should have length equal to input w")

  # Check inputs xy and d.
  if (!(length(xy) == p & length(d) == p))
    stop("Inputs xy and d should have length = ncol(X)")

  # Check inputs alpha0 and mu0.
  if (any(c(dim(alpha0),dim(mu0)) != c(p,K,p,K)))
    stop(paste("Inputs alpha0 and mu0 should be p x K matrices,",
               "with p = ncol(X) and K = length(w)"))
   
  # Check input Xr0.
  if (length(Xr0) != n)
    stop("length(Xr0) must be equal to nrow(X)")
   
  # Check input i.
  if (sum(i < 1 | i > p) > 0)
    stop("Input i contains invalid variable indices")

  # Initialize storage for the results.
  alpha <- t(alpha0)
  mu    <- t(mu0)
  Xr    <- c(Xr0)

  # Execute the C routine using the .Call interface and return the
  # updated variational parameters statistics in a list object. The
  # main reason for using the .Call interface is that there is less of
  # a constraint on the size of the input matrices. The only
  # components that change are alpha, mu and Xr. Note that I need to
  # subtract 1 from the indices because R vectors start at 1, and C
  # arrays start at 0. Also note that the alpha and mu matrices are
  # stored differently in the C implementation---variables correspond
  # to columns---so we need to first transpose these matrices.
  out <- .Call(C_varbvsmixupdate_Call,X = X,sigma = as.double(sigma),
               sa = as.double(sa),w = as.double(w),xy = as.double(xy),
               d = as.double(d),alpha = alpha,mu = mu,Xr = Xr,
               i = as.integer(i-1),eps = eps)
  return(list(alpha = t(alpha),mu = t(mu),Xr = Xr))
}
