% var1(X) produces the same result as var(X,1)', but does not require
% storage of any intermediate products of the same size as X.
function y = var1 (X)

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

  % Make sure X is in single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the sum of the variances.
  y = var1mex(X);
