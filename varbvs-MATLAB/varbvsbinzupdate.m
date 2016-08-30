% [alpha,mu,Xr] = varbvsbinzupdate(X,sa,logodds,stats,alpha0,mu0,Xr0,i) runs
% a single iteration of the coordinate ascent updates to maximize the
% variational lower bound for Bayesian variable selection in logistic
% regression, allowing for additional covariates. See varbvsbinupdate for
% more details.
function [alpha, mu, Xr] = ...
        varbvsbinzupdate (X, sa, logodds, stats, alpha0, mu0, Xr0, i)

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
      varbvsbinzupdatemex(X,double(sa),double(logodds),stats,double(alpha0),...
                          double(mu0),double(Xr0),double(i-1));
