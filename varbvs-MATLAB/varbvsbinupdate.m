% [alpha,mu,Xr] = varbvsbinupdate(X,sa,logodds,stats,alpha0,mu0,Xr0,i) runs
% a single iteration of the coordinate ascent updates to maximize the
% variational lower bound for Bayesian variable selection in logistic
% regression.
%
% Input X is an n x p matrix of observations of the variables (or features),
% where n is the number of samples, and p is the number of variables. Input
% y contains samples of the binary outcome; it is a vector of length n. X
% must be a single precision matrix.
%
% Input sa specifies the prior variance of the coefficients. Input logodds
% is the prior log-odds of inclusion for each variable. It must be a vector
% of length p. Note that a residual variance parameter (sigma) is not needed
% to model a binary outcome. See function updatestats in varbvsbin.m for
% more information about input 'stats'.
%
% Inputs alpha0, mu0 are the current parameters of the variational
% approximation; under the variational approximation, the ith regression
% coefficient is normal with probability alpha0(i), and mu0(i) is the mean
% of the coefficient given that it is included in the model. Input Xr0 must
% be Xr0 = X*(alpha0.*mu0).
%
% Input i specifies the order in which the coordinates are updated. It may
% be a vector of any length. Each entry of i must be an integer between 1
% and p.
%
% There are three outputs. Output vectors alpha and mu are the updated
% variational parameters, and Xr = X*(alpha.*mu). The computational
% complexity is O(n*length(i)).
function [alpha, mu, Xr] = ...
        varbvsbinupdate (X, sa, logodds, stats, alpha0, mu0, Xr0, i)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % X must be single precision.
  if ~isa(X,'single')
    error('Input X must be SINGLE');
  end

  % Check input sa.
  if ~isscalar(sa)
    error('Input sa must be a scalar');
  end

  % Check input logodds, alpha0 and mu0.
  if ~(length(logodds) == p & length(alpha0) == p & length(mu0) == p)
    error('logodds, alpha0 and mu0 must have length = size(X,2).');
  end

  % Check input Xr0.
  if length(Xr0) ~= n
    error('length(Xr0) must be equal to size(X,1)');
  end
   
  % Check input i.
  if sum(i < 1 | i > p)
    error('Input i contains invalid variable indices');
  end
   
  % Execute the C routine. I subtract 1 from the indices because MATLAB
  % arrays start at 1, but C arrays start at 0.
  [alpha mu Xr] = ...
      varbvsbinupdatemex(X,double(sa),double(logodds),stats,double(alpha0),...
                         double(mu0),double(Xr0),double(i-1));
