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
% Input X is an n x p matrix of observations about the variables (or
% features), where n is the number of samples, and p is the number of
% variables. Input y is the vector of observations about the outcome; it is
% a vector of length n.
%
% Inputs sigma, sa and logodds are additional model parameters; sigma and sa
% are scalars. Input sigma specifies the variance of the residual, and sa
% specifies the prior variance of the coefficients. logodds is the prior
% log-odds of inclusion for each variable. 
%
% Output logw is the variational estimate of the marginal log-likelihood
% given the hyperparameters. Outputs alpha, mu and s are the parameters of
% the variational approximation and, equivalently, variational estimates of
% posterior quantites: under the variational approximation, the ith
% regression coefficient is normal with probability alpha(i); mu(i) and s(i)
% are the mean and variance of the coefficient given that it is included in
% the model. 
%
% When update_sa = true, there is the additional option of computing the
% maximum a posteriori (MAP) estimate of the prior variance parameter (sa),
% in which sa is drawn from a scaled inverse chi-square distribution with
% scale sa0 and degrees of freedom n0.
function [logw, sigma, sa, alpha, mu, s] = ...
        varbvsnorm (X, y, sigma, sa, logodds, alpha, mu, tol, maxiter, ...
                    verbose, outer_iter, update_sigma, update_sa, n0, sa0)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);
  
  % (1) INITIAL STEPS
  % -----------------
  % Input X must be single precision.
  if ~isa(X,'single')
    X = single(X);
  end
  
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

    % Save the current variational parameters.
    alpha0 = alpha;
    mu0    = mu;

    % (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute the lower bound to the marginal log-likelihood.
    logw0 = int_linear(Xr,d,y,sigma,alpha,mu,s) ...
            + int_gamma(logodds,alpha) ...
            + int_klbeta(alpha,mu,s,sigma*sa);
    
    % (2b) UPDATE VARIATIONAL APPROXIMATION
    % -------------------------------------
    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      i = 1:p;
    else
      i = p:-1:1;
    end
    [alpha mu Xr] = varbvsnormupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,i);

    % (2c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute the lower bound to the marginal log-likelihood.
    logw = int_linear(Xr,d,y,sigma,alpha,mu,s) ...
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
    err = abs(alpha - alpha0);
    if verbose
      status = sprintf('%05d %05d %+13.6e %0.1e %06.1f %0.1e %0.1e',...
                       outer_iter,iter,logw,max(err),sum(alpha),sigma,sa);
      fprintf(status);
      fprintf(repmat('\b',1,length(status)));
    end
    if logw < logw0
      alpha = alpha0;
      mu    = mu0;
      logw  = logw0;
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
function I = int_linear (Xr, d, y, sigma, alpha, mu, s)
  n = length(y);
  I = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
      - d'*betavar(alpha,mu,s)/(2*sigma);
