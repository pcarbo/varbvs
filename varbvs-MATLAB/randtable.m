% randtable(p,n) returns n discrete random draws from the probability
% distribution with unnormalized contingency table p.
function x = randtable (p, n)

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

  % Get the number of random draws to generate.
  if ~exist('n')
    n = 1;
  end

  % Get the number of classes.
  m = numel(p);

  % Normalize the probability table.
  p = p(:)';
  p = p/sum(p);
  
  % Create the CDF.
  ub = cumsum(p);
  lb = [0 ub(1:m-1)];
  
  % Generate the discrete random deviates.
  x = zeros(n,1);
  for i = 1:n
    u    = repmat(rand,1,m);
    x(i) = sum((lb <= u & u < ub) .* (1:m));
  end
