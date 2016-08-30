% TO DO: Implement generalization of spike-and-slab with K normal
% components. For now, I assume no additional covariates Z. Note
% connection to "mvash" (special case in which we have individual-level
% data, and linear regression model).

% NOTES:
%
%   Options:
%     - tol
%     - maxiter
%     - verbose
%     - sigma
%     - q
%     - update_sigma
%     - update_q
%     - alpha
%     - mu
%
%   Is it important to have a penalty/regularization term for the
%   mixture weights (q)?
%
function fit = varbvsmix (X, Z, y, labels, sa, options)

% [logw,sigma,sa,alpha,mu,s] = varbvsnorm(X,y,sigma,sa,logodds,...)
% implements the fully-factorized variational approximation for Bayesian
% variable selection in linear regression. It finds the "best"
% fully-factorized variational approximation to the posterior distribution
% of the coefficients for a linear regression model of a continuous outcome
% (quantitiative trait), with spike and slab priors on the coefficients. By
% "best", we mean the approximating distribution that locally minimizes the
% K-L divergence between the approximating distribution and the exact
% posterior.
%
% Input X is an n x p matrix of variable (or feature) observations, where n
% is the number of samples, and p is the number of variables. Input y
% contains observations of the outcome; it is a vector of length n.
%
% Inputs sigma, sa and logodds are additional model parameters; sigma and sa
% are scalars. Input sigma specifies the variance of the residual, and sa
% specifies the prior variance of the coefficients (scaled by sigma). Input
% logodds is the prior log-odds of inclusion for each variable. Note that
% the prior log-odds here is defined with respect to the *natural*
% logarithm, whereas in function varbvs the prior log-odds is defined
% with respect to the base-10 logarithm, so a scaling factor of log(10)
% is needed to convert from the latter to the former.
%
% Output logw is the variational estimate of the marginal log-likelihood
% given the hyperparameters at each iteration of the co-ordinate ascent
% optimization procedure. Output err is the maximum difference between the
% approximate posterior probabilities (alpha) at successive iterations.
% Outputs alpha, mu and s are the parameters of the variational
% approximation and, equivalently, variational estimates of posterior
% quantites: under the variational approximation, the ith regression
% coefficient is normal with probability alpha(i); mu(i) and s(i) are the
% mean and variance of the coefficient given that it is included in the
% model.
%
% When update_sa = true, there is the additional option of computing the
% maximum a posteriori (MAP) estimate of the prior variance parameter (sa),
% in which sa is drawn from a scaled inverse chi-square distribution with
% scale sa0 and degrees of freedom n0.

  % Get the number of samples (n), variables (p) and mixture components (K).
  [n p] = size(X);
  K     = numel(sa);  
  
  % (1) CHECK INPUTS
  % ----------------
  % Input X must be single precision, and cannot be sparse.
  if issparse(X)
    error('Input X cannot be sparse');
  end
  if ~isa(X,'single')
    X = single(X);
  end

  % If input Z is not empty, it must be double precision, and must have as
  % many rows as X.
  if ~isempty(Z)
    if size(Z,1) ~= n
      error('Inputs X and Z do not match.');
    end
    Z = double(full(Z));
  end

  % Add intercept.
  Z    = [ones(n,1) Z];
  ncov = size(Z,2) - 1;

  % Input y must be a double-precision column vector with n elements.
  y = double(y(:));
  if length(y) ~= n
    error('Inputs X and y do not match.');
  end

  % The labels must be a cell array with p elements, or an empty array.
  if nargin < 4
    labels = [];
  end
  if isempty(labels)
    labels = cellfun(@num2str,num2cell(1:p)','UniformOutput',false);
  else
    labels = labels(:);
    if (~iscell(labels) | length(labels) ~= p)
      error('labels must be a cell array with numel(labels) = size(X,2)');
    end
  end
  
  % (2) PROCESS OPTIONS
  % -------------------
  % If the 'options' input argument is not specified, all the options are
  % set to the defaults.
  if nargin < 6
    options = [];
  end
  
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
  if ~isfinite(maxiter)
    error('options.maxiter must be a finite number');
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
    sigma = double(options.sigma);
  else
    sigma = var(y);
  end

  % OPTIONS.Q
  % Get the initial estimate of the mixture weights, if provided. By
  % default, the initial estimate is set so that all the weights are equal.
  if isfield(options,'q')
    q = double(options.q(:));
  else
    q = ones(K,1)/K;
  end

  % OPTIONS.UPDATE_SIGMA
  % Determine whether to update the residual variance parameter.
  if isfield(options,'update_sigma')
    update_sigma = options.update_sigma;
  else
    update_sigma = true;
  end

  % OPTIONS.UPDATE_Q
  % Determine whether to update the mixture weights.
  if isfield(options,'update_q')
    update_q = options.update_q;
  else
    update_q = true;
  end

  % OPTIONS.ALPHA
  % Set initial estimates of variational parameter alpha.
  if isfield(options,'alpha')
    alpha = double(options.alpha);
    if size(alpha,1) ~= p
      error('options.alpha must have one row for each variable (column of X)');
    end
    if size(alpha,2) ~= K
      error('options.alpha must have one column for each mixture component');
    end
  else
    alpha = rand(p,K);
    alpha = alpha ./ repmat(sum(alpha),p,1);
  end

  % OPTIONS.MU
  % Set initial estimates of variational parameter mu.
  if isfield(options,'mu')
    mu = double(options.mu);
    if size(mu,1) ~= p
      error('options.mu must have one row for each variable (column of X)');
    end
    if size(mu,2) == 1
      error('options.mu must have one column for each mixture component');
    end
  else
    mu = randn(p,K);
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
  SZX = (Z'*Z)\(Z'*X);
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
  
  % Compute a few useful quantities. Here I calculate X'*y as (y'*X)' to
  % avoid computing the transpose of X, since X may be large.
  xy = double(y'*X)';
  d  = diagsq(X);
  Xr = double(X*(alpha.*mu));
  
  % Calculate the variance of the coefficients.
  s = zeros(p,K);
  for i = 1:K
    s(:,i) = sigma*sa(i)./(sa(i)*d + 1);
  end
  
  % Initialize storage for outputs logw and err.
  logw = zeros(1,maxiter);
  err  = zeros(1,maxiter);
  
  % (4) MODEL FITTING - MAIN LOOP
  % -----------------------------
  % Repeat until convergence criterion is met, or until the maximum
  % number of iterations is reached.
  for iter = 1:maxiter

    % Save the current variational parameters and model parameters.
    alpha0 = alpha;
    mu0    = mu;
    s0     = s;
    sigma0 = sigma;
    q0     = q;

    % (4a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute the lower bound to the marginal log-likelihood.
    %
    % TO DO: Update this.
    %
    logw0 = p/2 - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
            - d'*betavarmix(alpha,mu,s)/(2*sigma);
    
    logw0 = - sum2(alpha.*log(alpha + eps))) + sum(alpha*log(q + eps)) ...
            I = (alpha'*log(s/sa) - alpha'*(s + mu.^2)/sa)/2 ...

            
    logw0 = int_linear(Xr,d,y,sigma,alpha,mu,s) ...
            + int_klbeta(alpha,mu,s,sigma*sa);
    
    % (4b) UPDATE VARIATIONAL APPROXIMATION
    % -------------------------------------
    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      i = 1:p;
    else
      i = p:-1:1;
    end
    % [alpha mu Xr] = varbvsnormupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,i);

    % (2c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute the lower bound to the marginal log-likelihood.
    logw(iter) = int_linear(Xr,d,y,sigma,alpha,mu,s) ...
                 + int_gamma(logodds,alpha) ...
                 + int_klbeta(alpha,mu,s,sigma*sa);
    
    % (2d) UPDATE RESIDUAL VARIANCE
    % -----------------------------
    % Compute the maximum likelihood estimate of sigma, if requested.
    % Note that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated.
    if update_sigma
      sigma = (norm(y - Xr)^2 + d'*betavar(alpha,mu,s) ...
               + alpha'*(s + mu.^2)/sa)/(n + sum(alpha));
      s = sa*sigma./(sa*d + 1);
    end
    
    % (2e) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    % -----------------------------------------------------
    % Compute the maximum a posteriori estimate of sa, if requested. Note
    % that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated. 
    if update_sa
      sa = (sa0*n0 + dot(alpha,s + mu.^2))/(n0 + sigma*sum(alpha));
      s  = sa*sigma./(sa*d + 1);
    end

    % (2f) CHECK CONVERGENCE
    % ----------------------
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum difference between the
    % posterior probabilities at two successive iterations is less than the
    % specified tolerance, or when the variational lower bound has
    % decreased. 
    err(iter) = max(abs(alpha - alpha0));
    if verbose
      if isempty(outer_iter)
        status = '';
      else
        status = sprintf('%05d ',outer_iter);
      end  
      status = [status sprintf('%05d %+13.6e %0.1e %06.1f %0.1e %0.1e',...
                               iter,logw(iter),err(iter),sum(alpha),sigma,sa)];
      fprintf(status);
      fprintf(repmat('\b',1,length(status)));
    end
    if logw(iter) < logw0
      logw(iter) = logw0;
      err(iter)  = 0;
      sigma      = sigma0;
      sa         = sa0;
      % alpha      = alpha0;
      mu         = mu0;
      s          = s0;
      break
    elseif err(iter) < tol
      break
    end
  end

  % Return the variational lower bound (logw) and "delta" in successive
  % iterates (err).
  logw = logw(1:iter);
  err  = err(1:iter);

% ----------------------------------------------------------------------
% Compute the lower bound to the marginal log-likelihood.
function I = compute_varlb (Xr, d, y, sigma, sa, q, alpha, mu, s)
  n = length(y);
  p = length(d);
  K = numel(sa);
  I = p/2 - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
          - d'*betavarmix(alpha,mu,s)/(2*sigma);
  for i = 1:K
    I = I - alpha(:,i)'*log(alpha(:,i) + eps) ...
          + sum(alpha(:,i)*log(q(i) + eps)) ...
          + alpha(:,i)'*log(s(:,i)/(sigma*sa(i)))/2 ...
          - alpha(:,i)'*(s(:,i) + mu(:,i).^2)/(sigma*sa(i))/2;
  end