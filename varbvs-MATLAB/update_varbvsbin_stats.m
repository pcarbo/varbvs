% This function is used in function 'varbvsbin' to compute quantities useful
% for the variational approximation to the logistic regression factors.
function stats = update_varbvsbin_stats (X, y, eta)

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

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute beta0 and yhat. See the journal paper for an explanation of
  % these two variables.
  beta0 = sum(y - 0.5)/sum(d);
  yhat  = y - 0.5 - beta0*d;

  % Calculate xy = X'*yhat as (yhat'*X)' and xd = X'*d as (d'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xd = double(d'*X)';

  % Compute the diagonal entries of X'*dhat*X. For a definition of dhat, see
  % the Bayesian Analysis journal paper.
  %
  % This is the less numerically stable version of this update:
  % 
  %   xdx = diagsq(X,d) - xd.^2/sum(d)
  % 
  dzr = d/sqrt(sum(d));
  xdx = diagsq(X,d) - double((dzr'*X).^2)';
 
  % Return the result.
  stats = struct('d',d,'yhat',yhat,'xy',xy,'xd',xd,'xdx',xdx);
