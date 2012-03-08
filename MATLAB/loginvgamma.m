% LOGINVGAMMA(X,A,B) returns the logarithm of the probability density
% function of the inverse gamma distribution at elements of X, with shape A
% and scale B.
function f = loginvgamma (x, a, b)
  f = a.*log(b) - gammaln(a) - (a+1).*log(x+eps) - b./(x+eps);
  