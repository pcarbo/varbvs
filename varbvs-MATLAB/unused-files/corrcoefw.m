% CORRCOEFW(X,Y,W) computes a Monte Carlo estimate of the correlation
% coefficient between random variables X and Y given samples of X and Y, and
% given normalized importance weights W. X, Y and W must be arrays of the
% same size.
function r = corrcoefw (x, y, w)
  r = covw(x,y,w) / sqrt(varw(x,w) * varw(y,w));
  