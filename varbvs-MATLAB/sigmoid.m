% SIGMOID(X) is the sigmoid of the elements of X. The sigmoid function is
% also known as the logistic link function. It is the inverse of LOGIT(X).
function y = sigmoid (x)
  y = 1./(1 + exp(-x));
