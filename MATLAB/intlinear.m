% INTLINEAR(XR,D,Y,SIGMA,A,MU,S) computes the expectation of the
% log-likelihood for the linear regression model of a quantitative trait,
% where the expectation is taken with respect to the fully-factorized
% variational approximation. Input Y is the vector of quantitative trait
% data. 
%
% Inputs A, MU and S are each vectors of the same length specifying the
% variational approximation; the Kth regression coefficient is normal with
% probability A(K), and zero with probability 1 - A(K). MU(K) and S(K) are
% the mean and variance of the normal mixture component. 
%
% XR must be equal to X*R, where R = ALPHA.*MU, and D = DIAG(X'*X).
function I = intlinear (Xr, d, y, sigma, alpha, mu, s)
  n = length(y);
  I = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
      - d'*betavar(alpha,mu,s)/(2*sigma);
