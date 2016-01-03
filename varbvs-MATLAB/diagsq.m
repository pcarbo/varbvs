% diagsq(X) is the same as diag(X'*X), but the computation is done more
% efficiently, and without having to store an intermediate matrix of the
% same size as X. diag(X,a) efficiently computes diag(X'*diag(a)*X). See the
% help for function diag for more details.
function y = diagsq (X, a)

  % If input a is not provided, set it to a vector of ones.
  [m n] = size(X);
  if ~exist('a')
    a = ones(m,1);
  end

  % Convert input a to a double column vector.
  a = double(a(:));

  % Convert X to single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the result using the efficient C routine.
  y = diagsqmex(X,a);
  