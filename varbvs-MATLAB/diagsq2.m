% diagsq2(X) is the same as diag(X'*A*X), but the computation is done
% more efficiently.
function y = diagsq2 (X, A)
  y = sum((X*A).*X,2);
    