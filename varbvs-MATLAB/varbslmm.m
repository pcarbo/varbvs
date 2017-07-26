%--------------------------------------------------------------------------
% varbslmm.m: Fit the "Bayesian sparse linear mixed model" (BSLMM) using
% variational approximation techniques.
%--------------------------------------------------------------------------
%
% TO DO: Add complete help here.
%
function fit = varbslmm (X, Z, y, logodds, sa, sb, options)

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
    
  % Get the number of samples (n), the number of variables (p), and the
  % number of hyperparameter settings (ns).
  [n p] = size(X);
  K     = numel(sa);

  % (1) CHECK INPUTS
  % ----------------
  % Input X should be single precision, and cannot be sparse.
  if issparse(X)
    error('Input X cannot be sparse');
  end
  if ~isa(X,'single')
    X = single(X);
  end

  % If input Z is not empty, it should be double precision, and should
  % have as many rows as X.
  if ~isempty(Z)
    if size(Z,1) ~= n
      error('Inputs X and Z do not match.');
    end
    Z = double(full(Z));
  end

  % Add the intercept.
  Z    = [ones(n,1) Z];
  ncov = size(Z,2) - 1;

  % Input y should be a double-precision column vector with n elements.
  y = double(y(:));
  if length(y) ~= n
    error('Inputs X and y do not match.');
  end

  % TO DO: Check that ns >= 2.

  % If the 'options' input argument is not specified, all the options are
  % set to the defaults.
  if nargin < 7
    options = [];
  end

  % (2) PROCESS OPTIONS
  % -------------------
  % OPTIONS.TOL
  % Set the convergence tolerance of the co-ordinate ascent updates.
  if isfield(options,'tol')
    tol = options.tol;
  else
    tol = 1e-4;
  end

  % OPTIONS.MAXITER
  % Set the maximum number of inner-loop iterations.
  if isfield(options,'maxiter')
    maxiter = options.maxiter;
  else
    maxiter = 1e4;
  end

  % OPTIONS.VERBOSE
  % Determine whether to output progress to the console.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end

  % OPTIONS.SIGMA
  % Get the initial estimate of the residual variance, if provided. By
  % default, the initial estimate is set to the sample variance of Y.
  if isfield(options,'sigma')
    sigma        = double(options.sigma);
    update_sigma = false;
  else
    sigma        = var(y);
    update_sigma = true;
  end
  if ~isscalar(sigma)
    error('options.sigma should be a scalar.')
  end

  % OPTIONS.UPDATE_SIGMA
  % Determine whether to update the residual variance parameter. Note
  % that the default setting is determined by whether options.sigma is
  % provided.
  if isfield(options,'update_sigma')
    update_sigma = options.update_sigma;
  end

  % OPTIONS.ALPHA
  % Set initial estimates of variational parameter alpha.
  initialize_params = true;
  if isfield(options,'alpha')
    alpha = double(options.alpha);
    initialize_params = false;  
    if size(alpha,1) ~= p
      error('options.alpha must have as many rows as X has columns');
    end
    if size(alpha,2) == 1
      alpha = repmat(alpha,1,ns);
    end
  else
    alpha = rand(p,ns);
    alpha = alpha ./ repmat(sum(alpha),p,1);
  end

  % OPTIONS.MU
  % Set initial estimates of variational parameter mu.
  if isfield(options,'mu')
    mu = double(options.mu);
    initialize_params = false;  
    if size(mu,1) ~= p
      error('options.mu must have as many rows as X has columns');
    end
    if size(mu,2) == 1
      mu = repmat(mu,1,ns);
    end
  else
    mu = randn(p,ns);
  end

  % (3) PREPROCESSING STEPS
  % -----------------------
  % Adjust the genotypes and phenotypes so that the linear effects of
  % the covariates are removed. This is equivalent to integrating out
  % the regression coefficients corresponding to the covariates with
  % respect to an improper, uniform prior; see Chipman, George and
  % McCulloch, "The Practical Implementation of Bayesian Model
  % Selection," 2001.
  %
  % Here I compute two quantities that are used here to remove linear
  % effects of the covariates (Z) on X and y, and later on (in function
  % "outerloop"), to efficiently compute estimates of the regression
  % coefficients for the covariates.
  SZy = (Z'*Z)\(Z'*y);
  SZX = double((Z'*Z)\(Z'*X));
  if ncov == 0
    X = X - repmat(mean(X),length(y),1);
    y = y - mean(y);
  else

    % This should give the same result as centering the columns of X and
    % subtracting the mean from y when we have only one covariate, the
    % intercept.
    y = y - Z*SZy;
    X = X - Z*SZX;
  end

  % Provide a brief summary of the analysis.
  if verbose
    % TO DO.
  end

  % (4) INITIALIZE STORAGE FOR THE OUTPUTS
  % --------------------------------------
  % TO DO.

  % (5) FIT BSLMM MODEL TO DATA
  % ---------------------------


  % (6) CREATE FINAL OUTPUT
  % -----------------------
  % TO DO.
