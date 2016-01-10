% var1(X) produces the same result as var(X,1)', but does not require
% storage of any intermediate products of the same size as X.
function y = var1 (X)
  
  % Make sure X is in single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
  % Compute the sum of the variances.
  y = var1mex(X);
