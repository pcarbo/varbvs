% diagsqt(X) is the same as diag(X*X'), but the computation is done more
% efficiently, and without having to store an intermediate result of the
% same size as X. diagsqt(X,a) efficiently computes diag(X*diag(a)*X').
function y = diagsqt (X, a)

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

  % If input a is not provided, set it to a vector of ones. Then convert
  % a to a double-precision column vector.
  [m n] = size(X);
  if ~exist('a')
    a = ones(n,1);
  end
  a = double(a(:));
  
  % Convert X to single precision, if necessary.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the result using the efficient C routine.
  y = diagsqtmex(X,a);
