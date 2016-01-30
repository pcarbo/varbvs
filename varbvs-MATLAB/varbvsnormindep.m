% TO DO: Add description of function (see varbvsnorm.m).
function [alpha, mu, s] = varbvsnormindep (X, y, sigma, sa, logodds)

  % Calculate the mean (mu) and variance (s) of the coefficients given that
  % they are included in the model, then calculate the posterior inclusion
  % probabilities (alpha), ignoring correlations between variables. Here I
  % calculate X'*y as (y'*X)' to avoid computing the transpose of X, since X
  % may be large.
  s     = sa*sigma./(sa*diagsq(X) + 1);
  mu    = s.*double(y'*X')/sigma;
  alpha = sigmoid(logodds + (log(s/(sa*sigma)) + mu.^2./s)/2);

