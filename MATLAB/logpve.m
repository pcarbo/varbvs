% LOGPVE(A,X) returns the logarithm of the density function for X, given
% that R = A*X/(A*X + 1) is uniform on [0,1]. Input A must be a positive
% scalar. This is useful for calculating the prior probability of the
% variance parameter X, when we have a uniform prior on R, the proportion of
% variance explained.
function f = logpve (a, x)
  r = a*x./(a*x+1);          % Proportion of variance explained.
  f = 2*log(r./x) - log(a);  % Log-density.
  