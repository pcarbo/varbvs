% *** DESCRIPTION OF FUNCTION GOES HERE ***
function I = intlogitz (stats, y)

  % Get some of the statistics.
  S = stats.S;
  d = stats.d;

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the approximate posterior distribution.
  I = sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 + logdet(S)/2 ...
      + qnorm(y - 0.5,S) u'*Z*(y - 0.5)/2 + yhat'*Xr - qnorm(Xr,U)^2/2 ...
      + a/2*(u'*Xr)^2 - d'*betavar(alpha,mu,s)/2;
