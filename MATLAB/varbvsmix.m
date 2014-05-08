% [LNZ,ALPHA,MU1,MU2,S1,S2,SIGMA] = VARBVSMIX(X,Y,SIGMA,SA1,SA2,LOGODDS)
% implements the fully-factorized variational approximation for the Bayesian
% variable selection mixture model with linear regression. It finds a "good"
% fully-factorized variational approximation to the posterior distribution
% of the coefficients in a linear regression model of a continuous outcome,
% with normal mixture priors on the coefficients. By "good", we mean the
% approximating distribution that locally minimizes the Kullback-Leibler
% divergence between the approximating distribution and the exact posterior.
%
% The required inputs are as follows. Input X is an N x P matrix of
% observations about the variables (or features), where N is the number of
% samples, and P is the number of variables. Y is the vector of observations
% about the outcome; it is a vector of length N. To account for an
% intercept, Y and X must be centered beforehand so that Y and each column
% of X has a mean of zero. Note that this function is not implemented quite
% as efficiently as VARBVS and VARBVSBIN to handle large matrices X.
%
% Inputs SA1, SA2 and LOGODDS are the fixed hyperparameters, in which SA and
% SA2 are scalars. SIGMA is a scalar specifying an initial estimate for the
% variance of the residual. SA1*SIGMA is the prior variance for coefficients
% drawn from the first mixture component, and SA2*SIGMA is the prior
% variance for coefficients drawn from the second mixture component. LOGODDS
% is the prior log-odds that a coefficient is drawn from the first mixture
% component. LOGODDS may either be a scalar, in which case all the variables
% have the same prior inclusion probability, or it may be a vector of length
% equal to P.
%
% There are seven outputs: output scalar LNZ is the variational estimate of
% the marginal log-likelihood given the hyperparameters; outputs ALPHA, MU1,
% MU2, S1 and S2 are the parameters of the variational approximation and,
% equivalently, variational estimates of posterior quantites. Outputs ALPHA,
% MU1, MU2, and S1 and S2 are all column vectors of length P. The final
% output is the maximum likelihood estimate of the residual variance, SIGMA.
%
% VARBVS(...,OPTIONS) overrides the default behaviour of the algorithm. Set
% OPTIONS.ALPHA, OPTIONS.MU1 and OPTIONS.MU2 to override the random
% initialization of variational parameters ALPHA, MU1 and MU2. And set
% OPTIONS.VERBOSE = FALSE to turn off reporting the algorithm's progress.
function [lnZ, alpha, mu1, mu2, s1, s2, sigma] = ...
    varbvsmix (X, y, sigma, sa1, sa2, logodds, options)

  % Convergence is reached when the maximum relative distance between
  % successive updates of the variational parameters is less than this
  % quantity.
  tolerance = 1e-4;  

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % CHECK INPUTS.
  % Note that I do not check the inputs as carefully as I do in VARBVS.
  y = y(:);

  % LOGODDS must be a double precision column vector of length P.
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  end
  logodds = logodds(:);

  % TAKE CARE OF OPTIONAL INPUTS.
  if ~exist('options')
    options = [];
  end

  % Set initial estimates of variational parameters.
  if isfield(options,'alpha')
    alpha = options.alpha(:);
  else
    alpha = rand(p,1);
    alpha = alpha / sum(alpha);
  end
  if isfield(options,'mu1')
    mu1 = options.mu1(:);
  else
    mu1 = sqrt(sigma*sa1) * randn(p,1);
  end
  if isfield(options,'mu2')
    mu2 = options.mu2(:);
  else
    mu2 = sqrt(sigma*sa2) * randn(p,1);
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
  xy = (y'*X)';
  d  = sum(X.^2)';
  Xr = X*betameanmix(alpha,mu1,mu2);
  
  % Calculate the variances of the coefficients.
  s1 = sa1*sigma./(sa1*d + 1);
  s2 = sa2*sigma./(sa2*d + 1);

  % MAIN LOOP.
  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  if verbose
    fprintf('       variational    max. incl max.      \n');
    fprintf('iter   lower bound  change vars E[b] sigma\n');
  end
  while true

    % Go to the next iteration.
    iter = iter + 1;

    % Save the current variational parameters and lower bound.
    alpha0  = alpha;
    mu10    = mu1;
    mu20    = mu2;
    lnZ0    = lnZ;
    params0 = [ alpha; alpha.*mu1; (1-alpha).*mu2 ];

    % UPDATE VARIATIONAL APPROXIMATION.
    % Run a forward or backward pass of the coordinate ascent updates.
    if isodd(iter)
      I = 1:p;
    else
      I = p:-1:1;
    end
    [alpha mu1 mu2 Xr] = varbvsmixupdate(X,sigma,sa1,sa2,logodds,xy,...
                                         d,alpha,mu1,mu2,Xr,I);

    % UPDATE RESIDUAL VARIANCE.
    % Compute the maximum likelihood estimate of the residual variance
    % parameter (SIGMA). Note that after updating the residual variance
    % parameter, we must also recalculate the variance of the regression
    % coefficients. 
    sigma = (norm(y - Xr)^2 + d'*betavarmix(alpha,mu1,mu2,s1,s2) ...
             + alpha'*(s1 + mu1.^2)/sa1 ...
             + (1-alpha)'*(s2 + mu2.^2)/sa2)/(n + p);
    s1 = sa1*sigma./(sa1*d + 1);
    s2 = sa2*sigma./(sa2*d + 1);

    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute the lower bound to the marginal log-likelihood.
    lnZ = intlinearmix(Xr,d,y,sigma,alpha,mu1,mu2,s1,s2) ...
          + intgamma(logodds,alpha) ...
          + intklbetamix(alpha,mu1,mu2,s1,s2,sigma*sa1,sigma*sa2);

    % CHECK CONVERGENCE.
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [ alpha; alpha.*mu1; (1-alpha).*mu2 ];
    I      = find(abs(params) > 1e-6);
    err    = relerr(params(I),params0(I));
    if verbose
      fprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f\n',iter,lnZ,max(err),...
              round(sum(alpha)),max(abs(alpha.*mu1)),sqrt(sigma));
    end
    if lnZ < lnZ0
      alpha = alpha0;
      mu1   = mu10;
      mu2   = mu20;
      lnZ   = lnZ0;
      break
    elseif max(err) < tolerance
      break
    end
  end