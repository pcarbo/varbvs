% slope(x) returns (sigmoid(x) - 1/2)./x, the slope of the conjugate to the
% log-sigmoid function at x, times 2. For details, see Bishop (2006), or the
% Bayesian Analysis paper. This is useful for working with the variational
% approximation for logistic regression.
function y = slope (x)
  y = (sigmoid(x) - 0.5)./(x + eps);
