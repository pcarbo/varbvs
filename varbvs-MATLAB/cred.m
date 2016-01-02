% [A,B] = CRED(X,W,X0) returns a 95% credible interval [A,B]. Precisely, we
% define the credible interval [A,B] to be the smallest interval containing
% X0 that contains 95% of the probability mass. (Note that the credible
% interval is not necessarily symmetric about X0. Other definitions of the
% credible interval are possible.) Input X are values of the random
% variable, and W contains the corresponding probabilities. These
% probabilities do not need to be normalized. X and W must be numeric arrays
% with the number of elements.
%
% [A,B] = CRED(X,W,X0,C) returns the credible interval for a number
% proportion different than 0.95. C must be between 0 and 1.
function [a, b] = cred (x, w, x0, c)

  % Set size of credible interval.
  if ~exist('c')
    c = 0.95;
  end

  % Convert the inputs X and W to column vectors.
  x = x(:);
  w = w(:);

  % Make sure the importance weights sum to 1.
  w = w / sum(w);

  % Get the number of points.
  n = length(x);

  % Sort the points in increasing order.
  [x I] = sort(x);
  w     = w(I);

  % Generate all possible intervals [a,b] from the set of points x.
  [a b] = ndgrid(1:n,1:n);
  I = find(a <= b);
  a = a(I);
  b = b(I);

  % Select only the intervals [a,b] that contain x0.
  I = find(x(a) <= x0 & x0 <= x(b));
  a = a(I);
  b = b(I);

  % Select only the intervals that contain 95% (or c%) of the mass.
  P = cumsum(w);
  I = find(P(b) - P(a) + w(a) >= c);
  a = a(I);
  b = b(I);

  % From the remaining intervals, choose the interval that has the
  % smallest width.
  [ans i] = min(x(b) - x(a));
  a = x(a(i));
  b = x(b(i));
