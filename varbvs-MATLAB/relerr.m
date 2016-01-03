% relerr(x,y) returns the absolute relative error between x and y. Note that
% relerr(x,y) = relerr(y,x), and we define relerr([],[]) = 0.
function y = relerr (x, y)
  if isempty(x) | isempty(y)
    y = 0;
  else
    y = abs(x-y)./(abs(x) + abs(y) + eps);
  end
  