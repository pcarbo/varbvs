% normalizelogweights(logw) takes as input an array of unnormalized
% log-probabilities logw and returns normalized probabilities such that the
% sum is equal to 1.
function w = normalizelogweights (logw)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %

  % Guard against underflow or overflow by adjusting the log-probabilities so
  % that the largest probability is 1.
  c = max(logw(:));
  w = exp(logw - c);

  % Normalize the probabilities.
  w = w / sum(w(:));  
  