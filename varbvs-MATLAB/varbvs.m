%--------------------------------------------------------------------------
% varbvs.m: One-sentence summary of function goes here.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Overview of function goes here.
%
% USAGE:
%    Summary of usage goes here.
%
% INPUT ARGUMENTS:
% Description of input arguments goes here.
%
% OUTPUT ARGUMENTS:
% Description of output arguments goes here.
%
% DETAILS:
%    Detailed description of function goes here.
%
% LICENSE: GPL v3
%
% DATE: December 28, 2015
%
% AUTHORS:
%    List contributors here.
%
% REFERENCES:
%    List of references goes here.
%
% SEE ALSO:
%    List related functions here.
%
% EXAMPLES:
%    Give some examples here.
%
function fit = varbvs (X, Z, y, family, options)

  % (1) CHECK INPUTS
  % ----------------
  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % Input X must be single precision.
  if ~isa(X,'single')
    X = single(X);
  end

  % If input Z is not empty, it must be double precision, and must have as
  % many rows as X.
  if ~isempty(Z)
    if size(Z,1) ~= n
      error('Inputs X and Z do not match.');
    end
    Z = double(Z);
  end

  % Input y must be a double-precision column vector with n elements.
  y = double(y(:));
  if length(y) ~= n
    error('Inputs X and y do not match.');
  end

  % By default, Y is a quantitative trait (normally distributed).
  if nargin < 4
    family = 'gaussian';
  end
  if isempty(family)
    family = 'gaussian';
  end
  if ~(family == 'gaussian' | family == 'binomial')
    error('family must be gaussian or binomial');
  end
  
  % (2) PROCESS OPTIONS
  % -------------------
  % If the 'options' input argument is not specified, all the options are
  % set to the defaults.
  if nargin < 5
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

  % OPTIONS.VERBOSE
  % Determine whether to output progress to the console.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end

  % OPTIONS.INTERCEPT
  % Determine whether to include an intercept in the model. If other
  % covariates are included, simply treat the intercept as another
  % covariate.
  if isfield(options,'intercept')
    intercept = options.intercept;
  else
    intercept = true;
  end
  if intercept & ~isempty(Z)
    Z = [ones(n,1) Z];
  end

  % OPTIONS.SIGMA
  % Get candidate settings for the variance of the residual, if provided.
  % Note that this option is not valid for a binary trait.
  if isfield(options,'sigma')
    sigma        = double(options.sigma(:)');
    update_sigma = false;
    if family == 'binomial'
      error('options.sigma is not valid with family = binomial')
    end
  else
    sigma        = var(y);
    update_sigma = true;
  end

  % OPTIONS.SA
  % Get candidate settings for the prior variance of the coefficients, if
  % provided.
  if isfield(options,'sa')
    sa        = double(options.sa(:)');
    update_sa = false;
  else
    sa        = 1;
    update_sa = true;
  end

  % OPTIONS.LOGODDS
  % Get candidate settings for the prior log-odds of inclusion. This may
  % either be specified as a vector, in which case this is the prior applied
  % uniformly to all variables, or it is a p x ns matrix, where p is the
  % number of variables and ns is the number of candidate hyperparameter
  % settings, in which case the prior log-odds is specified separately for
  % each variable. A default setting is only available if the number of
  % other hyperparameter settings is 1, in which case we select 20 candidate
  % settings for the prior log-odds, evenly spaced between log10(1/p) and
  % log10(0.5). If necessary, I convert the prior log-odds settings to an p
  % x ns matrix.
  if isfield(options,'logodds')
    logodds = double(options.theta(:)');
  elseif isscalar(sigma) & isscalar(sa)
    logodds = linspace(-log10(p),-0.3,20);
  else
    error('options.logodds must be specified')
  end
  if ~ismatrix(logodds) | size(logodds,1) ~= p
    logodds = repmat(logodds,p,1);
  end

  % Here is where I ensure that the numbers of candidate hyperparameter
  % settings are consistent.
  ns = max([numel(sigma) numel(sa) size(logodds,2)]);
  if isscalar(sigma)
    sigma = repmat(sigma,1,ns);
  end
  if isscalar(sa)
    sa = repmat(sa,1,ns);
  end
  if numel(sigma) ~= ns | numel(sa) ~= ns | size(logodds,2) ~= ns
    error('options.sigma, options.sa and options.logodds are inconsistent')
  end

  % OPTIONS.UPDATE_SIGMA
  % Determine whether to update the residual variance parameter. Note
  % that this option is not valid for a binary trait.
  if isfield(options,'update_sigma')
    update_sigma = options.update_sigma;
    if family == 'binomial'
      error('options.update_sigma is not valid with family = binomial')
    end
  end

  % OPTIONS.UPDATE_SA
  % Determine whether to update the prior variance of the regression
  % coefficients.
  if isfield(options,'update_sa')
    update_sa = options.update_sa;
  end

  % OPTIONS.SA0
  % Get the scale parameter for the scaled inverse chi-square prior.
  if isfield(options,'sa0')
    sa0 = options.sa0;
  else
    sa0 = 0;
  end

  % OPTIONS.N0
  % Get the number of degrees of freedom for the scaled inverse chi-square
  % prior.
  if isfield(options,'n0')
    n0 = options.n0;
  else
    n0 = 0;
  end

  % OPTIONS.ALPHA
  % Set initial estimates of variational parameter alpha.
  if isfield(options,'alpha')
    alpha = double(options.alpha);
    initialize_params = false;  
    if size(alpha,1) ~= p
      error('options.alpha must have as many rows as X has columns')
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
    mu = double(options.mu(:));
    initialize_params = false;  
    if size(mu,1) ~= p
      error('options.mu must have as many rows as X has columns')
    end
    if size(mu,2) == 1
      mu = repmat(mu,1,ns);
    end
  else
    mu = randn(p,ns);
  end

  % OPTIONS.ETA
  % Set initial estimates of variational parameter eta. Note this is only
  % relevant for logistic regression.
  if isfield(options,'eta')
    eta = double(options.eta);
    initialize_params = false;  
    if family ~= 'binomial'
      error('options.eta is only valid for family = binomial');
    end
    if size(eta,1) ~= n
      error('options.eta must have as many rows as X')
    end
    if size(eta,2) == 1
      eta = repmat(eta,1,ns);
    end
  elseif family == 'binomial'
    eta = ones(n,ns);
  end

  % OPTIONS.INITIALIZE_PARAMS
  % Determine whether to find a good initialization for the variational
  % parameters.
  if isfield(options,'initialize_params')
    initialize_params = options.initialize_params;
  end
  
  % TO DO: Allow specification of summary statistics from "fixed"
  % variational estimates for an external set of variables.
  clear options

  % (3) PREPROCESSING STEPS
  % -----------------------
  % TO DO.
  
  % (4) INITIALIZE STORAGE FOR THE OUTPUTS
  % --------------------------------------
  % Initialize storage for the variational estimate of the marginal
  % log-likelihood for each hyperparameter setting (logw), and the variances
  % of the regression coefficients (s).
  logw = zeros(1,ns);
  s    = zeros(p,ns);
  
  % (5) FIT BAYESIAN VARIABLE SELECTION MODEL TO DATA
  % -------------------------------------------------
  if ns == 1

    % TO DO: Implement special case when there is only one hyperparameter
    % setting.
    
  else
      
    % If a good initialization isn't already provided, find a good
    % initialization for the variational parameters. Repeat for each
    % candidate setting of the hyperparameters.
    if initialize_params
      fprintf('Finding best initialization for %d combinations of ',ns);
      fprintf('hyperparameters.\n');
      for i = 1:ns

        % Find a set of parameters that locally minimize the Kullback-Leibler
        % divergence between the approximating distribution and the exact
        % posterior.
        if family == 'gaussian'
        end
      end
    
      % Choose an initialization common to all the runs of the coordinate
      % ascent algorithm. This is chosen from the hyperparameters with
      % the highest variational estimate of the marginal likelihood.
      [ans i] = max(logw);
      alpha   = repmat(alpha(:,i),1,ns);
      mu      = repmat(mu(:,i),1,ns);
      eta     = repmat(eta(:,i),1,ns);
    end
    
    % Compute a variational approximation to the posterior distribution
    % for each candidate setting of the hyperparameters.
    fprintf('Computing marginal likelihood for %d combinations of ',ns);
    fprintf('hyperparameters.\n');
    for i = 1:ns

      % Find a set of parameters that locally minimize the Kullback-Leibler
      % divergence between the approximating distribution and the exact
      % posterior.
      % TO DO.
    end
  end

  % 5. CREATE FINAL OUTPUT
  % ----------------------
  if family == 'gaussian'
    fit = struct('logw',logw,'sigma',sigma,'sa',sa,'logodds',logodds,...
                 'alpha',alpha,'mu',mu,'s',s);
  elseif family == 'binomial'
    fit = struct('logw',logw,'sigma',sigma,'sa',sa,'logodds',logodds,...
                 'alpha',alpha,'mu',mu,'s',s,'eta',eta);
  end
  