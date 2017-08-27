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
# TO DO: Add brief description of function here.
varbslmm <- function (X, Z, y, sigma, logodds, h, d, update.sigma,
                      tol = 1e-4, maxiter = 1e4, verbose = TRUE) {

  # (1) CHECK INPUTS
  # ----------------
  # Check input X.
  if (!(is.matrix(X) & is.double(X) & sum(is.na(X)) == 0))
    stop("Input X must be a double-precision matrix with no missing values.")
    
  # Add column names to X if they are not provided.
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X",1:p)

  # Check input Z.
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) | sum(is.na(Z)) > 0)
      stop("Input Z must be a numeric matrix with no missing values.")
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    storage.mode(Z) <- "double"
  }

  # Add intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)

  # Check input y.
  if (!is.numeric(y) | sum(is.na(y)) > 0)
    stop("Input y must be a numeric vector with no missing values.")
  if (length(y) != n)
    stop("Inputs X and y do not match")
  y <- c(as.double(y))

  
}
