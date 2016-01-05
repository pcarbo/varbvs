% VARW(X,W) computes a Monte Carlo estimate of the variance of X given
% samples of X and normalized importance weights W. X and W must be arrays
% of the same size.
function v = varw (x, w)
  v = covw(x,x,w);
