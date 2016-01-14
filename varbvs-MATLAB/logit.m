% logit(x) returns the logit of the elements of X. It is the inverse of
% sigmoid(x).
function y = logit (x)
  y = log((x + eps)./((1 - x) + eps));
