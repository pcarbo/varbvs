% LOGIT(X) is the logit of the elements of X. This the inverse of
% SIGMOID(X).
function y = logit (x)
  y = log((x + eps)./((1 - x) + eps));
