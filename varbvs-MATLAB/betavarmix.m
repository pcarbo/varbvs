% betavarmix(p,mu,s) returns variances of variables drawn from mixtures of
% normals. Each of the inputs is a n x k matrix, where n is the number of
% variables and k is the number of mixture components. Specifically,
% variable i is drawn from a mixture in which the jth mixture component is
% the univariate normal with mean mu(i,j) and variance s(i,j).
%
% Note that the following two lines should return the same result when k=2
% and the first component is the "spike" density with zero mean and
% variance.
%
%   y1 = betavar(p,mu,s)
%   y2 = betavarmix([1-p p],[zeros(n,1) mu],[zeros(n,1) s])
%
function y = betavarmix (p, mu, s)

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
  y = sum(p.*(s + mu.^2),2) - sum(p.*mu,2).^2;
