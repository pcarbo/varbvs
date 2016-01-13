% logpexp(x) returns log(1 + exp(x)). The computation is performed in a
% numerically stable manner. For large entries of x, log(1 + exp(x)) is
% effectively the same as x.
function y = logpexp (x)
  y    = x;
  i    = find(x < 16);
  y(i) = log(1 + exp(x(i)));
