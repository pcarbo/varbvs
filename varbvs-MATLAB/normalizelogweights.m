% normalizelogweights(logw) takes as input an array of unnormalized
% log-probabilities logw and returns normalized probabilities such that the
% sum is equal to 1.
function w = normalizelogweights (logw)

  % Guard against underflow or overflow by adjusting the log-probabilities so
  % that the largest probability is 1.
  c = max(logw(:));
  w = exp(logw - c);

  % Normalize the probabilities.
  w = w / sum(w(:));  
  