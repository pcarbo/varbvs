% logdet(X) safely calculates of logarithm of the determinant of X when X is
% positive definite.
function y = logdet (X)
  y = 2*sum(log(diag(chol(X))));
