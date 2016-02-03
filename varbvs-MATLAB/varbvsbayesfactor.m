
% Compute the marginal log-likelihood under the null hypothesis using
% importance sampling.
c     = max(logw0(:));
logZ0 = c + log(mean(exp(logw0(:) - c)));

% Compute the marginal log-likelihood under the enrichment hypothesis using
% importance sampling.
c     = max(logw1(:));
logZ1 = c + log(mean(exp(logw1(:) - c)));

  % Get the numerical estimate of the Bayes factor.
  BF = exp(logZ1 - logZ0);
