% LOGBETA(X,A,B) returns the logarithm of the probability density function
% of the Beta distribution at elements of X, with prior sample sizes A and B.
function f = logbeta (x, a, b)
  f = (a-1).*log(x+eps) + (b-1).*log(1-x+eps) - betaln(a,b);
  