% diagsqt(X) is the same as diag(X*X'), but the computation is done more
% efficiently, and without having to store an intermediate result of the
% same size as X. diagsqt(X,a) efficiently computes diag(X*diag(a)*X').
function y = diagsqt (X, a)

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
