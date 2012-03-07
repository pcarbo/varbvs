% RELERR(X,Y) returns the absolute relative error between X and Y. They must
% be array of the same dimension. Note that RELERR(X,Y) = RELERR(Y,X). We
% define RELERR([],[]) = 0.
function y = relerr (x1, x2)
  if isempty(x1) & isempty(x2)
    y = 0;
  else
    x = (abs(x1) + abs(x2));
    y = abs(x1 - x2) ./ (x + eps);
  end
  