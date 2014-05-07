% *** FIX THESE COMMENTS ***
% 
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
function [lnZ, alpha, mu1, mu2, s1, s2] = ...
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
  Xr = X*(alpha.*mu1 + (1-alpha).*mu2);
  
  % Calculate the variances of the coefficients.
  s1 = sa1*sigma./(sa1*d + 1);
  s2 = sa2*sigma./(sa2*d + 1);

  % MAIN LOOP.
  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  fprintf('       variational    max. incl max.      \n');
  fprintf('iter   lower bound  change vars E[b] sigma\n');
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
    for i = 1:p
  
      % Update the variational estimates of the posterior means.
      r      = alpha(i)*mu1(i) + (1-alpha(i))*mu2(i);
      mu1(i) = s1(i)/sigma * (xy(i) + d(i)*r - dot(X(:,i),Xr));
      mu2(i) = s2(i)/sigma * (xy(i) + d(i)*r - dot(X(:,i),Xr));
  
      % Update the variational estimate of the posterior inclusion
      % probability.
      SSR1     = mu1(i)^2/s1(i);
      SSR2     = mu2(i)^2/s2(i);
      alpha(i) = sigmoid(logodds + (log(s1(i)*sa2/(s2(i)*sa1)) + SSR1-SSR2)/2);
  
      % Update Xr = X*r.
      rnew = alpha(i)*mu1(i) + (1-alpha(i))*mu2(i);
      Xr   = Xr + (rnew - r)*X(:,i);
    end

    % UPDATE RESIDUAL VARIANCE.
    sigma = (norm(y - Xr)^2 + d'*betavarmix(alpha,mu1,mu2,s1,s2) ...
             + alpha'*(s1 + mu1.^2)/sa1 ...
             + (1-alpha)'*(s2 + mu2.^2)/sa2)/(n + p);
    s1 = sa1*sigma./(sa1*d + 1);
    s2 = sa2*sigma./(sa2*d + 1);

    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute the lower bound to the marginal log-likelihood.
    lnZ = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
          - d'*betavarmix(alpha,mu1,mu2,s1,s2)/(2*sigma) ...
          + intgamma(logodds,alpha) ...
	  + intklbeta(alpha,mu1,s1,sigma*sa1) ...
          + intklbeta(1-alpha,mu2,s2,sigma*sa2) ...
          + alpha'*log(alpha + eps) + (1 - alpha)'*log(1 - alpha + eps);

    % CHECK CONVERGENCE.
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [ alpha; alpha.*mu1; (1-alpha).*mu2 ];
    I      = find(abs(params) > 1e-6);
    err    = relerr(params(I),params0(I));
    fprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f\n',iter,lnZ,max(err),...
            round(sum(alpha)),max(abs(alpha.*mu1)),sqrt(sigma));
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