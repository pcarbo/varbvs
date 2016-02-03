% TO DO: Describe this function here.
function BF = varbvsbayesfactor (logw0, logw1)

  % Compute the marginal log-likelihood under the null hypothesis using
  % importance sampling.
  c     = max(logw0(:));
  logz0 = c + log(mean(exp(logw0(:) - c)));

  % Compute the marginal log-likelihood under the alternative hypothesis
  % using importance sampling.
  c     = max(logw1(:));
  logz1 = c + log(mean(exp(logw1(:) - c)));

  % Compute the numerical estimate of the Bayes factor.
  BF = exp(logz1 - logz0);
