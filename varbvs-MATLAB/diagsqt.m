% DIAGSQT(X) is the same as DIAG(X*X'), but the computation is done more
% efficiently, and without having to store an intermediate result of the
% same size as X. 
%
% DIAGSQT(X,A) efficiently computes DIAG(X*DIAG(A)*X').
function y = diagsqt (X, a)

  % If input A is not provided, set it to a vector of ones.
  [m n] = size(X);
  if ~exist('a')
    a = ones(n,1);
  end
  if length(a) ~= n
    error('Inputs X and A do not match');
  end

  % Convert A to a double column vector.
  a = double(a(:));
  
  % Convert X to single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the result using the efficient C routine.
  y = diagsqtmatlab(X,a);

  