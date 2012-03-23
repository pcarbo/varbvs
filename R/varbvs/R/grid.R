grid3d <- function (x, y, z) {
  # Does the same thing as ndgrid(x,y,z) in MATLAB.
  
  # Get the number of entries in each of the inputs.
  nx <- length(x)
  ny <- length(y)
  nz <- length(z)

  # Initialize the 3-d outputs.
  X <- array(dim=c(nx,ny,nz))
  Y <- X
  Z <- X

  # Set the entries of the 3-D arrays.
  for (i in 1:nx)
    for (j in 1:ny)
      for (k in 1:nz) {
        X[i,j,k] <- x[i]
        Y[i,j,k] <- y[j]
        Z[i,j,k] <- z[k]
      }
  
  return(list(X = X,Y = Y,Z = Z))
}
