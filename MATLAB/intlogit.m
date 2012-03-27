% INTLOGIT(Y,STATS,ALPHA,MU,S,XR,ETA) computes an integral that appears in
% the variational lower bound of the marginal log-likelihood for the
% logistic regression model. This integral is an approximation to the
% expectation of the logistic regression log-likelihood taken with respect
% to the variational approximation. 
%
% Input arguments Y and STATS specify the variational approximation to the
% likelihood. Y is the column vector of observations about the binary
% outcome. STATS is the STRUCT output from UPDATESTATS. Inputs ALPHA, MU and
% S are the parameters of the approximating distribution; see VARBVSBIN for
% details. ETA is the vector of free parameters specifying the variational
% lower bound to the logistic regression factors. Input XR must be equal to
% XR = X*R, where R is the posterior mean of the coefficients. Note that
% under the fully-factorized variational approximation, R = ALPHA.*MU; see
% VARBVSBIN for details. ETA and XR are column vectors with the same length
% as Y.
function I = intlogit (y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  d    = stats.d;
  u    = stats.u;
  U    = diag(sparse(u));

  % Get the variance of the intercept given the other coefficients.
  a = 1/sum(u);

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the variational approximation.
  I = sum(logsigmoid(eta)) + eta'*(u.*eta - 1)/2 + log(a)/2 ...
      + a*sum(y - 0.5)^2/2 + yhat'*Xr - qnorm(Xr,U)^2/2 ...
      + a/2*(u'*Xr)^2 - d'*betavar(alpha,mu,s)/2;
