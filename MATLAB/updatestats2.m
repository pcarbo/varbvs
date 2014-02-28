function varargout = updatestats (X, Z, y, eta)

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute the posterior covariance of u (the regression coefficients
  % corresponding to the Z variables) given beta (the regression
  % coefficients corresponding to the X variables).
  