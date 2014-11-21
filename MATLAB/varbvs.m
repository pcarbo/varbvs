% [LNZ,ALPHA,MU,S] = VARBVS(X,Y,SIGMA,SA,LOGODDS) implements the
% fully-factorized variational approximation for Bayesian variable selection
% in linear regression. It finds the "best" fully-factorized variational
% approximation to the posterior distribution of the coefficients in a
% linear regression model of a continuous outcome (quantitiative trait),
% with spike and slab priors on the coefficients. By "best", we mean the
% approximating distribution that locally minimizes the Kullback-Leibler
% divergence between the approximating distribution and the exact posterior.
%
% The required inputs are as follows. Input X is an N x P matrix of
% observations about the variables (or features), where N is the number of
% samples, and P is the number of variables. Y is the vector of observations
% about the outcome; it is a vector of length N. To account for an
% intercept, Y and X must be centered beforehand so that Y and each column
% of X has a mean of zero. 
%
% Note that this routine is implemented with the assumption that the data X
% is single floating-point precision (type HELP SINGLE), as opposed to the
% MATLAB default of double precision. This is useful for large data sets,
% because single precision requires half of the number of bits as double
% floating-point precision. If X is provided in another numerical
% representation, it is immediately converted to SINGLE.
%
% Inputs SIGMA, SA and LOGODDS are the hyperparameters. SIGMA and SA are
% scalars. SIGMA specifies the variance of the residual, and SA*SIGMA is the
% prior variance of the coefficients. LOGODDS is the prior log-odds of
% inclusion for each variable. It is equal to LOGODDS = LOG(Q./(1-Q)), where
% Q is the prior probability that each variable is included in the linear
% model of Y (i.e. the prior probability that its coefficient is not zero).
% LOGODDS may either be a scalar, in which case all the variables have the
% same prior inclusion probability, or it may be a vector of length P.
%
% There are four outputs. Output scalar LNZ is the variational estimate of
% the marginal log-likelihood given the hyperparameters SIGMA, SA and
% LOGODDS.
%
% Outputs ALPHA, MU and S are the parameters of the variational
% approximation and, equivalently, variational estimates of posterior
% quantites: under the variational approximation, the ith regression
% coefficient is normal with probability ALPHA(i); MU(i) and S(i) are the
% mean and variance of the coefficient given that it is included in the
% model. Outputs ALPHA, MU and S are all column vectors of length P.
%
% VARBVS(...,OPTIONS) overrides the default behaviour of the algorithm. Set
% OPTIONS.ALPHA and OPTIONS.MU to override the random initialization of
% variational parameters ALPHA and MU. Set OPTIONS.VERBOSE = FALSE to turn
% off reporting the algorithm's progress. Set OPTIONS.UPDATE_VARS to specify
% the indices of the variables to update, so that the variational estimates
% for the other variables remain fixed. And set OPTIONS.UPDATE_SIGMA = TRUE
% to compute the maximum likelihood estimate of the residual variance
% (SIGMA), in which case input SIGMA acts as the initial estimate of this
% parameter. When this option is set to TRUE, the maximum likelihood
% estimate is returned in the fifth output.
function [lnZ, alpha, mu, s, sigma, sa] = ...
        varbvs (X, y, sigma, sa, logodds, options)

  % Convergence is reached when the maximum relative distance between
  % successive updates of the variational parameters is less than this
  % quantity.
  tolerance = 1e-4;  

  % CHECK INPUTS.
  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % X must be single precision.
  if ~isa(X,'single')
    X = single(X);
  end

  % Y must be a double precision column vector of length N.
  y = double(y(:));
  if length(y) ~= n
    error('Data X and Y do not match');
  end

  % SIGMA and SA must be double scalars.
  sigma = double(sigma);
  sa    = double(sa); 
  if ~(isscalar(sigma) & isscalar(sa))
    error('Inputs SIGMA and SA must be scalars');
  end

  % LOGODDS must be a double precision column vector of length P.
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  end
  logodds = double(logodds(:));
  if length(logodds) ~= p
    error('LOGODDS must be a scalar or a vector of length P');
  end

  % TAKE CARE OF OPTIONAL INPUTS.
  if ~exist('options')
    options = [];
  end

  % Set initial estimates of variational parameters.
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

  % Determine whether to update the variational estimates for only a
  % subset of the variables.
  if isfield(options,'update_vars')
    update_vars = options.update_vars(:)';
  else
    update_vars = 1:p;
  end
  p = length(update_vars);

  % Determine whether to update the residual variance parameter.
  if isfield(options,'update_sigma')
    update_sigma = options.update_sigma;
  else
    update_sigma = false;
  end

  % Determine whether to update the parameter specifying the prior variance of
  % the additive effects.
  if isfield(options,'update_sa')
    update_sa = options.update_sa;
  else
    update_sa = false;
  end

  % Determine whether to display the algorithm's progress.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end
  clear options
  
  % INITIAL STEPS.
  % Compute a few useful quantities. Here I calculate X'*Y as (Y'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(y'*X)';
  d  = diagsq(X);
  Xr = double(X*(alpha.*mu));
  
  % Calculate the variance of the coefficients.
  s = sa*sigma./(sa*d + 1);
  
  % MAIN LOOP.
  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  if verbose
    fprintf('       variational    max. incl max.           \n');
    fprintf('iter   lower bound  change vars E[b] sigma   sd\n');
  end
  while true

    % Go to the next iteration.
    iter = iter + 1;
    
    % Save the current variational parameters and lower bound.
    alpha0  = alpha;
    mu0     = mu;
    lnZ0    = lnZ;
    params0 = [ alpha; alpha .* mu ];

    % UPDATE VARIATIONAL APPROXIMATION.
    % Run a forward or backward pass of the coordinate ascent updates.
    if isodd(iter)
      I = update_vars(1:p);
    else
      I = update_vars(p:-1:1);
    end
    [alpha mu Xr] = varbvsupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,I);
    
    % UPDATE RESIDUAL VARIANCE.
    % Compute the maximum likelihood estimate of the residual variance
    % parameter (SIGMA), if requested. Note that we must also recalculate
    % the variance of the regression coefficients.
    if update_sigma
      sigma = (norm(y - Xr)^2 + d'*betavar(alpha,mu,s) ...
               + alpha'*(s + mu.^2)/sa)/(n + sum(alpha));
      s = sa*sigma./(sa*d + 1);
    end

    % UPDATE PRIOR VARIANCE OF ADDITIVE EFFECTS.
    % Compute the maximum likelihood estimate of the prior variance of the
    % additive effects (SA), if requested. Note that we must also
    % recalculate the variance of the regression coefficients.
    if update_sa
      sa = dot(alpha,s + mu.^2)/(sigma*sum(alpha));
      s  = sa*sigma./(sa*d + 1);
    end

    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute the lower bound to the marginal log-likelihood.
    lnZ = intlinear(Xr,d,y,sigma,alpha,mu,s) ...
	  + intgamma(logodds,alpha) ...
	  + intklbeta(alpha,mu,s,sigma*sa);
    
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
      fprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f %0.2f\n',iter,lnZ,max(err),...
	      round(sum(alpha)),max(abs(alpha.*mu)),sqrt(sigma),sqrt(sa));
    end
    if lnZ < lnZ0
      alpha = alpha0;
      mu    = mu0;
      lnZ   = lnZ0;
      break
    elseif max(err) < tolerance
      break
    end
  end
