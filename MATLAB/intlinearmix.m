% INTLINEAR(XR,D,Y,SIGMA,ALPHA,MU,S) computes an integral that appears in
% the variational lower bound of the marginal log-likelihood for the mixture
% model. This integral is the expectation of the linear regression
% log-likelihood taken with respect to the variational approximation to the
% mixture model.
function I = intlinearmix (Xr, d, y, sigma, alpha, mu1, mu2, s1, s2)
  n = length(y);
  I = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
      - d'*betavarmix(alpha,mu1,mu2,s1,s2)/(2*sigma);
