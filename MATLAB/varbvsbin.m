% [LNZ,ALPHA,MU,S,ETA] = VARBVSBIN(X,Y,SA,LOGODDS) implements the
% fully-factorized variational approximation for Bayesian variable selection
% in logistic regression. It finds the "best" fully-factorized variational
% approximation to the posterior distribution of the coefficients in a
% logistic regression model of a binary outcome or trait (e.g. disease
% status in a case-control study), with spike and slab priors on the
% coefficients. By "best", we mean the approximating distribution that
% locally minimizes the Kullback-Leibler divergence between the
% approximating distribution and the exact posterior.
%
% The required inputs are as follows. Input X is an N x P matrix of
% observations about the variables (or features), where N is the number of
% samples, and P is the number of variables. Y is the vector of observations
% about the binary trait; it is a vector of length N. Unlike function
% VARBVS, Y and X must *not* be centered. Instead, we will account for the
% intercept as we update the variational approximation.
%
% Note that this routine is implemented with the assumption that the data X
% is single floating-point precision (type HELP SINGLE), as opposed to the
% MATLAB default of double precision. This is useful for large data sets,
% because single precision requires half of the number of bits as double
% floating-point precision. If X is provided in another numerical
% representation, it is immediately converted to SINGLE.
%
% Inputs SA and LOGODDS are the hyperparameters. Scalar SA is the prior
% variance of the coefficients. LOGODDS is the prior log-odds of inclusion
% for each variable. It is equal to LOGODDS = LOG(Q./(1-Q)), where Q is the
% prior probability that each variable is included in the linear model of Y
% (i.e. the prior probability that its coefficient is not zero). LOGODDS may
% either be a scalar, in which case all the variables have the same prior
% inclusion probability, or it may be a vector of length P. Note that the
% residual variance parameter (SIGMA) is not needed to model a binary trait.
%
% There are five outputs. Output scalar LNZ is the variational estimate of
% the marginal log-likelihood given the hyperparameters SA and LOGODDS.
%
% Outputs ALPHA, MU and S are the parameters of the variational
% approximation and, equivalently, variational estimates of posterior
% quantites: under the variational approximation, the ith regression
% coefficient is normal with probability ALPHA(i); MU(i) and S(i) are the
% mean and variance of the coefficient given that it is included in the
% model. Outputs ALPHA, MU and S are all column vectors of length P.
%
% ETA is the vector of free parameters that specifies the variational
% approximation to the likelihood factors in the logistic regression. ETA is
% a vector of length N.
%
% VARBVS(...,OPTIONS) overrides the default behaviour of the algorithm. Set
% OPTIONS.ALPHA and OPTIONS.MU to override the random initialization of
% variational parameters ALPHA and MU. Set OPTIONS.ETA to override the
% default initialization of the free parameters specifying the variational
% lower bound on the logistic regression factors. Set OPTIONS.FIXED_ETA =
% TRUE to prevent ETA from being updated. And set OPTIONS.VERBOSE = FALSE to
% turn off reporting the algorithm's progress.
function [lnZ, alpha, mu, s, eta] = varbvsbin (X, y, sa, logodds, options)

  % Convergence is reached when the maximum relative distance between
  % successive updates of the variational parameters is less than this
  % quantity.
  tolerance = 1e-4;  
  
  % CHECK INPUTS.
  % Get the number of samples and variables.
  [n p] = size(X);

  % Make sure the genotypes are in single precision.
  if ~isa(X,'single')
    X = single(X);
  end

  % Y must be a double precision column vector of length N.
  y = double(y(:));
  if length(y) ~= n
    error('Data X and Y do not match');
  end

  % SA must be a double scalar.
  sa = double(sa); 
  if ~isscalar(sa)
    error('Input SA must be scalars');
  end

  % LOGODDS must be a double precision column vector of length P.
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  end
  logodds = double(logodds(:));
  if length(logodds) ~= p
    error('LOGODDS must be a scalar or a vector of the same length as Y');
  end

  % TAKE CARE OF OPTIONS.
  if ~exist('options')
    options = [];
  end

  % Set initial estimates of sufficient statistics.
  if isfield(options,'alpha') 
    alpha = double(options.alpha(:));
  else
    alpha = rand(p,1);
    alpha = alpha / sum(alpha);
  end
  if isfield(options,'mu')
    mu = double(options.mu(:));
  else
    mu = randn(p,1);
  end
  if length(alpha) ~= p || length(mu) ~= p
    error('OPTIONS.ALPHA and OPTIONS.MU must be vectors of length P');
  end

  % Initialize the free parameters specifying the variational approximation
  % to the logistic regression factors.
  if isfield(options,'eta')
    eta = options.eta;
  else
    eta = ones(n,1);
  end

  % Determine whether to update the variational approximation to the
  % logistic regression.
  if isfield(options,'fixed_eta')
    if ~isfield(options.eta)
      error('OPTIONS.FIXED_ETA = TRUE requires input OPTIONS.ETA');
    end
    fixed_eta = options.fixed_eta;
  else
    fixed_eta = false;
  end
  
  % Determine whether to display the algorithm's progress.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end
  clear options
  
  % INITIAL STEPS.
  % Compute a few useful quantities.
  Xr    = double(X*(alpha.*mu));
  stats = updatestats(X,y,eta);

  % MAIN LOOP.
  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  if verbose
    fprintf('       variational    max. incl max.\n');
    fprintf('iter   lower bound  change vars E[b]\n');
  end
  while true

    % Go to the next iteration.
    iter = iter + 1;
    
    % Save the current variational parameters and lower bound.
    alpha0  = alpha;
    mu0     = mu;
    lnZ0    = lnZ;
    eta0    = eta;
    params0 = [ alpha; alpha .* mu ];

    % UPDATE VARIATIONAL APPROXIMATION.
    % Run a forward or backward pass of the coordinate ascent updates.
    if isodd(iter)
      I = 1:p;
    else
      I = p:-1:1;
    end
    [alpha mu Xr] = varbvsbinupdate(X,sa,logodds,stats,alpha,mu,Xr,I);

    % Recalculate the posterior variance of the coefficients.
    s = sa./(sa*stats.d + 1);  

    % UPDATE ETA.
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    if ~fixed_eta
      eta   = updateeta(X,y,betavar(alpha,mu,s),Xr,stats.u);
      stats = updatestats(X,y,eta);
    end

    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute variational lower bound to marginal log-likelihood.
    lnZ = intlogit(y,stats,alpha,mu,s,Xr,eta) ...
	  + intgamma(logodds,alpha) ...
	  + intklbeta(alpha,mu,s,sa);
    
    % CHECK CONVERGENCE.
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [ alpha; alpha .* mu ];
    I      = find(abs(params) > 1e-6);
    err    = relerr(params(I),params0(I));
    if verbose
     fprintf('%4d %+13.6e %0.1e %4d %0.2f\n',iter,lnZ,max(err),...
	      round(sum(alpha)),max(abs(alpha.*mu)));
    end
    if lnZ < lnZ0
      alpha = alpha0;
      mu    = mu0;
      eta   = eta0;
      lnZ   = lnZ0;
      break
    elseif max(err) < tolerance
      break
    end
  end
