% int_klbeta(alpha,mu,s,sa) computes an integral that appears in the
% variational lower bound of the marginal log-likelihood. This integral is
% the negative K-L divergence between the approximating distribution and the
% prior of the coefficients. Note that this sa is not the same as the sa
% used as an input to varbvsnorm.
function I = int_klbeta (alpha, mu, s, sa)

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
  I = (sum(alpha) + alpha'*log(s/sa) - alpha'*(s + mu.^2)/sa)/2 ...
      - alpha'*log(alpha + eps) - (1 - alpha)'*log(1 - alpha + eps);
  