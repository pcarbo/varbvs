% Use this instead of log(sigmoid(x)) to avoid loss of numerical precision.
function y = logsigmoid (x)
  y = -logpexp(-x);
  