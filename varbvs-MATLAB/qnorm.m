% QNORM(X,A) returns the quadratic norm of vector X with respect to positive
% definite matrix A. For a definition of the quadratic norm, see p. 635 of
% Convex Optimization (2004) by Boyd & Vandenberghe.
function y = qnorm (x, A)
  x = x(:);
  y = sqrt(x'*A*x);
  