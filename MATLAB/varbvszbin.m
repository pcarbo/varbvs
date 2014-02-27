% *** DESCRIPTION OF FUNCTION GOES HERE ***
function [lnZ, alpha, mu, s, eta] = varbvszbin (X, Z, y, sa, logodds, ...
                                                options)

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
    if ~isfield(options,'eta')
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
    % TO DO.

    % Recalculate the posterior variance of the coefficients.
    % TO DO.

    % UPDATE ETA.
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    if ~fixed_eta
      % TO DO.
    end
    
    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute variational lower bound to marginal log-likelihood.
    % TO DO.

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
