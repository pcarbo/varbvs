% TO DO: Update this function description.
%
% [logw,ALPHA,MU,S,SIGMA,SA] = varbvsnorm(X,Y,SIGMA,SA,LOGODDS,...) implements the
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
% for the other variables remain fixed. Set OPTIONS.UPDATE_SIGMA = TRUE
% to compute the maximum likelihood estimate of the residual variance
% (SIGMA), in which case input SIGMA acts as the initial estimate of this
% parameter. When this option is set to TRUE, the maximum likelihood
% estimate is returned as one of the outputs. Set OPTIONS.UPDATE_SA =
% TRUE to compute the maximum likelihood estimate of the prior variance
% of the regression coefficients (SA), in which case input SA acts as the
% initial estimate of this parameter. When this option is set to TRUE,
% the maximum likelihood estimate is returned as one of the outputs.
%
% When OPTIONS.UPDATE_SA = TRUE, there is the additional option of computing
% the maximum a posteriori estimate of the prior variance parameter (SA),
% in which SA is drawn from a scaled inverse chi-square distribution with
% scale OPTIONS.SA0 and degrees of freedom OPTIONS.N0.
function [logw, sigma, sa, alpha, mu, s] = ...
        varbvsnorm (X, y, sigma, sa, logodds, tol, maxiter, verbose, ...
                    outer_iter, update_sigma, update_sa, n0, sa0)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);
  
  % (1) INITIAL STEPS
  % -----------------
  % Compute a few useful quantities. Here I calculate X'*y as (y'*X)' to
  % avoid computing the transpose of X, since X may be large.
  xy = double(y'*X)';
  d  = diagsq(X);
  Xr = double(X*(alpha.*mu));
  
  % Calculate the variance of the coefficients.
  s = sa*sigma./(sa*d + 1);
  
  % (2) MAIN LOOP
  % -------------
  % Repeat until convergence criterion is met, or until the maximum
  % number of iterations is reached.
  logw = -Inf;
  for iter = 1:maxiter

    % Save the current variational parameters and lower bound.
    alpha0  = alpha;
    mu0     = mu;
    logw0   = logw;
    params0 = [ alpha; alpha .* mu ];

    % (2a) UPDATE VARIATIONAL APPROXIMATION
    % -------------------------------------
    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      i = 1:p;
    else
      i = p:-1:1;
    end
    [alpha mu Xr] = varbvsupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,i);
    
    % (2b) UPDATE RESIDUAL VARIANCE
    % -----------------------------
    % Compute the maximum likelihood estimate of sigma, if requested.
    % Note that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated.
    if update_sigma
      sigma = (norm(y - Xr)^2 + d'*betavar(alpha,mu,s) ...
               + alpha'*(s + mu.^2)/sa)/(n + sum(alpha));
      s = sa*sigma./(sa*d + 1);
    end

    % (2c) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    % -----------------------------------------------------
    % Compute the maximum a posteriori estimate of sa, if requested. Note
    % that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated. 
    if update_sa
      sa = (sa0*n0 + dot(alpha,s + mu.^2))/(n0 + sigma*sum(alpha));
      s  = sa*sigma./(sa*d + 1);
    end

    % (2d) COMPUTE VARIATIONAL LOWER BOUND
    % ------------------------------------
    % Compute the lower bound to the marginal log-likelihood.
    logw = intlinear(Xr,d,y,sigma,alpha,mu,s) ...
           + intgamma(logodds,alpha) ...
           + intklbeta(alpha,mu,s,sigma*sa);
    
    % (2e) CHECK CONVERGENCE
    % ----------------------
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small. If the variational bound decreases,
    % stop.
    params = [alpha; alpha.*mu];
    i      = find(abs(params) > 1e-6);
    err    = relerr(params(i),params0(i));
    if verbose
      % TO DO: Fix this.
      fprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f %0.2f\n',iter,logw,...
              max(err),round(sum(alpha)),max(abs(alpha.*mu)),sigma,sa);
    end
    if lnZ < lnZ0
      alpha = alpha0;
      mu    = mu0;
      lnZ   = lnZ0;
      break
    elseif max(err) < tol
      break
    end
  end

% ----------------------------------------------------------------------
% intlinear(Xr,d,y,sigma,alpha,mu,s) computes an integral that appears in
% the variational lower bound of the marginal log-likelihood. This integral
% is the expectation of the linear regression log-likelihood taken with
% respect to the variational approximation.
function y = intlinear (Xr, d, y, sigma, alpha, mu, s)
  n = length(y);
  y = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
      - d'*betavar(alpha,mu,s)/(2*sigma);
