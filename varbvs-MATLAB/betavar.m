% betavar(p,mu,s) returns the variance of X, in which X is drawn from the
% normal distribution with probability p, and X is 0 with probability
% 1-p. Inputs mu and s specify the mean and variance of the normal
% density. Inputs p, mu and s must be arrays of the same dimension. This
% function is useful for calculating the variance of the coefficients under
% the fully-factorized variational approximation.
function y = betavar (p, mu, s)

  % Note that this is the same as 
  % 
  %    v = p*(s + mu^2) - (p*mu)^2.
  %
  y = p.*(s + (1 - p).*mu.^2);
