% int_gamma(logodds,alpha) computes an integral that appears in the
% variational lower bound of the marginal log-likelihood. This integral is
% the expectation on the prior inclusion probabilities taken with respect to
% the variational approximation.
function I = int_gamma (logodds, alpha)

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

  % This is the same as 
  %
  %    sum(alpha.*log(q) + (1-alpha).*log(1-q)).
  %  
  I = sum((alpha-1).*logodds + logsigmoid(logodds));
