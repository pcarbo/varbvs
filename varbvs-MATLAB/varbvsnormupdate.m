% [alpha,mu,Xr] = varbvsnormupdate(X,sigma,sa,logodds,xy,d,alpha0,mu0,Xr0,i)
% runs a single iteration of the coordinate ascent updates to maximize the
% variational lower bound for Bayesian variable selection in linear
% regression. 
%
% Input X is an n x p matrix of variable (or feature) observations, where n
% is the number of samples, and p is the number of variables. Input xy =
% X'*y, where y is the vector of samples of the continuous outcome. X must
% be a single precision matrix.
%
% Inputs sigma, sa and logodds specify other model parameters. sigma and sa
% are scalars. sigma specifies the variance of the residual, and sa is the
% prior variance of the regression coefficients (scaled by sigma). Input
% logodds is the prior log-odds of inclusion for each variable. It must be a
% vector of length p.
%
% Inputs alpha0, mu0 are the current parameters of the variational
% approximation. Under the variational approximation, the ith regression
% coefficient is included in the model with probability alpha0(i), and
% mu0(i) is the mean of the coefficient given that it is included in the
% model. Inputs Xr0 and d must be Xr0 = X*(alpha0.*mu0) and d = diag(X'*X).
%
% Input i specifies the order in which the coordinates are updated. It may
% be a vector of any length. Each entry of i must be an integer between 1
% and p.
%
% There are 3 outputs. Output vectors alpha and mu are the updated
% variational parameters, and Xr = X*(alpha.*mu). The computational
% complexity is O(n*length(i)).
function [alpha, mu, Xr] = ...
        varbvsnormupdate (X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, i)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % X must be single precision.
  if ~isa(X,'single')
    error('Input X must be SINGLE');
  end

  % Check inputs sigma and sa.
  if ~isscalar(sigma) | ~isscalar(sa)
    error('Inputs sigma and sa must be scalars');
  end
   
  % Check input logodds, xy, d, alpha0 and mu0.
  if ~(length(logodds) == p & length(xy) == p & length(d) == p & ...
       length(alpha0) == p & length(mu0) == p)
    error('logodds, xy, d, alpha0 and mu0 must have length = size(X,2).');
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
      varbvsnormupdatemex(X,double(sigma),double(sa),double(logodds),...
                          double(xy),double(d),double(alpha0),double(mu0),...
                          double(Xr0),double(i-1));