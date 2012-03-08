% DIAGSQ(X) is the same as DIAG(X'*X), but the computation is done more
% efficiently, and without having to store an intermediate matrix of the
% same size as X. 
%
% DIAGSQ(X,A) efficiently computes DIAG(X'*DIAG(A)*X).
function y = diagsq (X, a)

  if ~exist('a')
    n = size(X,1);
    a = ones(n,1);
  end

  % Convert X to single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the result using the C++ routine.
  y = diagsqfast(X,a);
  