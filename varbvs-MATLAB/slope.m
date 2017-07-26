% slope(x) returns (sigmoid(x) - 1/2)./x, the slope of the conjugate to the
% log-sigmoid function at x, times 2. For details, see Bishop (2006), or the
% Bayesian Analysis paper. This is useful for working with the variational
% approximation for logistic regression.
function y = slope (x)

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
  y = (sigmoid(x) - 0.5)./(x + eps);
