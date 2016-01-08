% [a,b] = cred(x,w,x0,c) returns a c% credible interval [a,b], in which c is
% a number between 0 and 1. Precisely, we define the credible interval [a,b]
% to be the smallest interval containing x0 that contains c% of the
% probability mass. (Note that the credible interval is not necessarily
% symmetric about x0. Other definitions of the credible interval are
% possible.) By default, c = 0.95.
%
% Input x is the vector of random variable assignments, and w contains the
% corresponding probabilities. (These probabilities need not be normalized.)
% Inputs x and w must be numeric arrays with the same number of elements.
function [a, b] = cred (x, w, x0, c)

  % Set the size of the credible interval.
  if nargin < 4
    c = 0.95;
  end

  % Convert the inputs x and w to column vectors.
  x = x(:);
  w = w(:);

  % Make sure the importance weights sum to 1.
  w = w/sum(w);

  % Get the number of points.
  n = length(x);

  % Sort the points in increasing order.
  [x i] = sort(x);
  w     = w(i);

  % Generate all possible intervals [a,b] from the set of points x.
  [a b] = ndgrid(1:n,1:n);
  i     = find(a <= b);
  a     = a(i);
  b     = b(i);

  % Select only the intervals [a,b] that contain x0.
  i = find(x(a) <= x0 & x0 <= x(b));
  a = a(i);
  b = b(i);

  % Select only the intervals that c% of the mass.
  P = cumsum(w);
  i = find(P(b) - P(a) + w(a) >= c);
  a = a(i);
  b = b(i);

  % From the remaining intervals, choose the interval that has the
  % smallest width.
  [ans i] = min(x(b) - x(a));
  a = x(a(i));
  b = x(b(i));
