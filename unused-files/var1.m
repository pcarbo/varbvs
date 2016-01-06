% VAR1(X) is the same as VAR(X,1)'.
function y = var1 (X)
  
  % Make sure X is in single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the sum of the variances.
  y = var1matlab(X);