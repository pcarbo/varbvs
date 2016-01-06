% COVW(X,Y,W) computes a Monte Carlo estimate of the covariance between
% random variables X and Y given samples of X and Y, and given normalized
% importance weights W. X, Y and W must be arrays of the same size.
function c = covw (x, y, w)

  % Convert the inputs to vectors.
  x = x(:);
  y = y(:);
  w = w(:);

  % Compute the mean of X and Y.
  muX = dot(w,x);
  muY = dot(w,y);

  % Compute the covariance between X and Y.
  c = dot(w,x.*y) - muX * muY;