% LOGIT10(X) is the logit of the elements of X in base 10 instead of base
% e. It is the inverse of SIGMOID10(X).
function y = logit10 (x)
  y = log10((x + eps)./((1 - x) + eps));
