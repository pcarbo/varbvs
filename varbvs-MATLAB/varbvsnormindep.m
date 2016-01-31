% This function computes the mean (mu) and variance (s) of the coefficients
% given that they are included in the linear regression model, then it
% computes the posterior inclusion probabilities (alpha), ignoring
% correlations between variables. This function is used in varbvsindep.m.
function [alpha, mu, s] = varbvsnormindep (X, y, sigma, sa, logodds)
  s     = sa*sigma./(sa*diagsq(X) + 1);
  mu    = s.*double(y'*X)'/sigma;
  alpha = sigmoid(logodds + (log(s/(sa*sigma)) + mu.^2./s)/2);

