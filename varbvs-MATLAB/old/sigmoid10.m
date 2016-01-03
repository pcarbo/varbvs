% SIGMOID10(X) is the sigmoid of the elements of X in base 10 instead of
% base e. It is the inverse of LOGIT10(X).
function y = sigmoid10 (x)
  y = 1./(1 + 10.^(-x));
