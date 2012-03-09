% *** FIX COMMENTS ***
% 
% INTLOGIT(Y,STATS,ALPHA,MU,S,XR,ETA) returns the variational lower bound on
% the expectation of the log-likelihood for the logistic regression model of
% a binary trait, where the expectation is taken with respect to the
% fully-factorized variational approximation. Input arguments Y and STATS
% specify the likelihood; Y is the vector of quantitative trait data, and
% for details on STATS see the help for function UPDATESTATS.
%
% Inputs ALPHA, MU and S are the parameters of the variational
% approximation; the Kth regression coefficient is normal with probability
% A(K), and zero with probability 1 - A(K). MU(K) and S(K) are the mean and
% variance of the normal mixture component. ALPHA, MU and S are vectors with
% one entry for every SNP. ETA is the vector of parameters of the
% variational approximation to the nonlinear terms in the logistic
% regression. And input XR should be equal to X*R, where R = ALPHA.*MU. ETA
% and XR are vectors of the same length as Y.
function I = intlogit (y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  d    = stats.d;
  u    = stats.u;
  U    = spdiag(u);

  % Get the variance of the intercept given the additive effects.
  a = 1/sum(u);

  % Compute the variational lower bound on the expectation of the
  % log-likelihood with respect to the variational approximation.
  I = sum(logsigmoid(eta)) + eta'*(u.*eta - 1)/2 + log(a)/2 ...
      + a*sum(y - 0.5)^2/2 + yhat'*Xr - qnorm(Xr,U)^2/2 ...
      + a/2*(u'*Xr)^2 - d'*betavar(alpha,mu,s)/2;
