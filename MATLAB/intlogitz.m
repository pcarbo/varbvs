% INTLOGITZ(Z,Y,STATS,ALPHA,MU,S,XR,ETA) is the same as INTLOGIT, except
% that it allows for an arbitrary set of covariates Z. STATS is the STRUCT
% output from UPDATESTATS2.
function I = intlogitz (Z, y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  xdx  = stats.xdx;
  S    = stats.S;
  d    = stats.d;
  D    = diag(sparse(d));

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the approximate posterior distribution.
  I = sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 + logdet(S)/2 ...
      + qnorm(Z'*(y - 0.5),S)^2/2 + yhat'*Xr - qnorm(Xr,D)^2/2 ...
      + qnorm(Z'*D*Xr,S)^2/2 - xdx'*betavar(alpha,mu,s)/2;
