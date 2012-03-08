% NORMALIZELOGWEIGHTS(LOGW) takes as input an array of unnormalized
% log-importance weights LOGW and returns normalized importance weights such
% that the sum of the normalized importance weights is equal to one.
function w = normalizelogweights (logw)

  % We guard against underflow or overflow by adjusting the log-importance
  % weights so that the largest importance weight is one.
  c = max(logw(:));
  w = exp(logw - c);

  % Normalize the importance weights.
  w = w / sum(w(:));  
  