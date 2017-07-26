% betavar(p,mu,s) returns the variance of X, in which X is drawn from the
% normal distribution with probability p, and X is 0 with probability
% 1-p. Inputs mu and s specify the mean and variance of the normal
% density. Inputs p, mu and s must be arrays of the same dimension. This
% function is useful for calculating the variance of the coefficients under
% the fully-factorized variational approximation.
function y = betavar (p, mu, s)

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
      
  % Note that this is the same as 
  % 
  %    v = p*(s + mu^2) - (p*mu)^2.
  %
  y = p.*(s + (1 - p).*mu.^2);
