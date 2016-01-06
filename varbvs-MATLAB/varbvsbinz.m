% [logw,sa,alpha,mu,s,eta] = varbvsbinz(X,Z,y,sa,logodds,...) implements the
% fully-factorized variational approximation for Bayesian variable selection
% in logistic regression, allowing for covariates. It is the same as
% varbvsbin, except that it allows for an additional set of covariates that
% are not subject to the same "spike-and-slab" priors as the other
% variables. The covariate data Z are specified as an n x m matrix, where n
% is the number of samples, and m is the number of covariates. This function
% is equivalent to varbvsbin when only one covariate is specified, the
% intercept, and Z = ones(n,1).
function [logw, sa, alpha, mu, s, eta] = ...
        varbvsbinz (X, Z, y, sa, logodds, alpha, mu, eta, tol, maxiter, ...
                   verbose, outer_iter, update_sa, optimize_eta, n0, sa0)
  
  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % (1) INITIAL STEPS
  % -----------------
  % Input X must be single precision.
  if ~isa(X,'single')
    X = single(X);
  end

  % Compute a few useful quantities.
  Xr    = double(X*(alpha.*mu));
  stats = update_stats(X,Z,y,eta);
  s     = sa./(sa*stats.xdx + 1);

  % (2) MAIN LOOP
  % -------------
  % Repeat until convergence criterion is met, or until the maximum
  % number of iterations is reached.
  logw = -Inf;
  for iter = 1:maxiter
    
    % Save the current variational parameters.
    alpha0  = alpha;
    mu0     = mu;
    eta0    = eta;
    params0 = [alpha; alpha.*mu];

    % (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    % --------------------------------------------
    logw0 = int_logit(Z,y,stats,alpha,mu,s,Xr,eta) ...
            + int_gamma(logodds,alpha) ...
            + int_klbeta(alpha,mu,s,sa);

    % (2b) UPDATE VARIATIONAL APPROXIMATION
    % -------------------------------------
    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      i = 1:p;
    else
      i = p:-1:1;
    end
    [alpha mu Xr] = varbvsbinzupdate(X,sa,logodds,stats,alpha,mu,Xr,i);

    % (2c) UPDATE ETA
    % ---------------
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    if optimize_eta
      eta   = update_eta(X,Z,y,betavar(alpha,mu,s),Xr,stats.d);
      stats = update_stats(X,Z,y,eta);
      s     = sa./(sa*stats.xdx + 1);
    end
   
    % (2d) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute variational lower bound to marginal log-likelihood.
    logw = int_logit(Z,y,stats,alpha,mu,s,Xr,eta) ...
           + int_gamma(logodds,alpha) ...
           + int_klbeta(alpha,mu,s,sa);

    % (2e) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    % -----------------------------------------------------
    % Compute the maximum a posteriori estimate of sa, if requested. Note
    % that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated. 
    if update_sa
      sa = (sa0*n0 + dot(alpha,s + mu.^2))/(n0 + sigma*sum(alpha));
      s  = sa./(sa*stats.xdx + 1);
    end

    % (2f) CHECK CONVERGENCE
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
      % TO DO: FIX THIS.
      fprintf('%4d %+13.6e %0.1e %4d %0.2f %0.2f\n',iter,logw,max(err),...
              round(sum(alpha)),max(abs(alpha.*mu)),sqrt(sa));
    end
    if logw < logw0
      alpha = alpha0;
      mu    = mu0;
      eta   = eta0;
      logw  = logw0;
      break
    elseif max(err) < tolerance
      break
    end
  end

% ----------------------------------------------------------------------
% diagprod(A,B) efficiently computes diagprod(A*B').
function y = diagprod (A, B)
  y = double(sum(A.*B,2));
  
% ----------------------------------------------------------------------
% update_stats(X,Z,y,eta) returns useful quantities for updating the
% variational approximation to the logistic regression factors, allowing for
% covariates.
function stats = update_stats (X, Z, y, eta)

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute the posterior covariance of u (coefficients for Z) given beta
  % (coefficients for X).
  D = diag(sparse(d));
  S = inv(Z'*D*Z);
  
  % Compute matrix dzr = D*Z*R', where R is an upper triangular matrix such
  % that R'*R = S.
  R   = chol(S);
  dzr = D*(Z*R');

  % Compute yhatT. 
  yhat = y - 0.5 - dzr*R*(Z'*(y - 0.5));

  % Here, I calculate xy = X'*yhat as (yhat'*X)' and xd = X'*d as (d'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xd = double(d'*X)';

  % Compute the diagonal entries of X'*Dhat*X. For a definition of Dhat,
  % see the Bayesian Analysis journal paper.
  xdx = diagsq(X,d) - diagsq(dzr'*X);

  % Return the result.
  stats = struct('S',S,'d',d,'yhat',yhat,'xy',xy,'xd',xd,'xdx',xdx,'dzr',dzr);

% ----------------------------------------------------------------------
% updateeta(X,Z,y,v,Xr,d) returns the M-step update for the parameters
% specifying the variational lower bound to the logistic regression factors,
% allowing for additional covariates.
function eta = update_eta (X, Z, y, v, Xr, d)

  % Compute MUZ, the posterior mean of u, the regression coefficients
  % corresponding to the covariates. Here, S is the posterior covariance of
  % u given beta.
  D   = diag(sparse(d));
  S   = inv(Z'*D*Z);
  muz = S*Z'*(y - 0.5 - d.*Xr);

  % Calculate the covariance between the coefficients u and beta.
  V = diag(sparse(v));
  C = -double((S*Z'*D)*X)*V;

  % This is the M-step update for the free parameters.
  U   = double((S*Z'*D)*X)';
  eta = sqrt((Z*muz + Xr).^2 + diagsq2(Z,S) + diagsq2(Z,U'*V*U) ...
             + diagsqt(X,v) + 2*diagprod(X*C',Z));

% -----------------------------------------------------------------------
function I = int_logit (Z, y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  xdx  = stats.xdx;
  S    = stats.S;
  d    = stats.d;
  D    = diag(sparse(d));

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the approximate posterior distribution.
  I = sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 + logdet(S)/2 ...
      + qnorm(Z'*(y - 0.5),S)^2/2 + yhat'*Xr - qnorm(Xr,D)^2/2 ...
      + qnorm(Z'*D*Xr,S)^2/2 - xdx'*betavar(alpha,mu,s)/2;
