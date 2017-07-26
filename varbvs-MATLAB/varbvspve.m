% varbvspve(X,fit,nr) samples nr posterior estimates of the proportion of
% variance in Y explained by the Bayesian variable selection model fitted
% using a variational approximation. This function is only valid for the
% linear regression model (family = 'gaussian') with an intercept.
function pve = varbvspve (X, fit, nr)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %

  % Get the number of variables.
  p = size(X,2);
    
  % Input X must be single precision, and cannot be sparse.
  if issparse(X)
    error('Input X cannot be sparse');
  end
  if ~isa(X,'single')
    X = single(X);
  end
  if (numel(fit.labels) ~= p)
    error('Inputs X and fit are not compatible');
  end
  
  % Take care of the optional inputs.
  if nargin < 3
    nr = 1000;
  end

  % Initialize storage for posterior estimates of the proportion of variance
  % explained.
  pve = zeros(nr,1);

  % For each sample, compute the proportion of variance explained.
  for i = 1:nr

    % Draw a hyperparameter setting from the posterior distribution.
    j = randtable(fit.w);
    
    % Sample the region coefficients.
    b = fit.mu(:,j) + sqrt(fit.s(:,j)) .* randn(p,1);
    b = b .* (rand(p,1) < fit.alpha(:,j));

    % Compute the proportion of variance explained.
    sz     = var(X*b,1);
    pve(i) = sz/(sz + fit.sigma(j));
  end  
