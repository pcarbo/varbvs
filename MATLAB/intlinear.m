% INTLINEAR(XR,D,Y,SIGMA,ALPHA,MU,S) computes an integral that appears in
% the variational lower bound of the marginal log-likelihood. This integral
% is the expectation of the linear regression log-likelihood taken with
% respect to the variational approximation. Inputs XR and D must be XR =
% X*(ALPHA.*MU) and D = DIAG(X'*X). For a description of the remaining
% inputs, see the help for function VARBVS.
function I = intlinear (Xr, d, y, sigma, alpha, mu, s)
  n = length(y);
  I = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
      - d'*betavar(alpha,mu,s)/(2*sigma);
