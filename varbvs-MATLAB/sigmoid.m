% sigmoid(x) returns the sigmoid of the elements of x. The sigmoid function
% is also known as the logistic link function. It is the inverse of
% logit(x).
function y = sigmoid (x)
  y = 1./(1 + exp(-x));
