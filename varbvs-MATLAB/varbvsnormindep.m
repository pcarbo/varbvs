% This function computes the mean (mu) and variance (s) of the coefficients
% given that they are included in the linear regression model, then it
% computes the posterior inclusion probabilities (alpha), ignoring
% correlations between variables. This function is used in varbvsindep.m.
function [alpha, mu, s] = varbvsnormindep (X, y, sigma, sa, logodds)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %
  s     = sa*sigma./(sa*diagsq(X) + 1);
  mu    = s.*double(y'*X)'/sigma;
  alpha = sigmoid(logodds + (log(s/(sa*sigma)) + mu.^2./s)/2);

