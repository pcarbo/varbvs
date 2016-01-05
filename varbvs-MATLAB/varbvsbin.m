% TO DO: Update description of this function.
% 
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
function [lnZ, alpha, mu, s, eta] = ...
        varbvsbin (X, y, sa, logodds, alpha, mu, eta, tol, maxiter, ...
                   verbose, outer_iter, update_sa, update_eta, n0, sa0)

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
  stats = update_stats(X,y,eta);

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
    eta0    = eta;
    params0 = [alpha; alpha.*mu];

    % (2a) UPDATE VARIATIONAL APPROXIMATION
    % -------------------------------------
    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      i = 1:p;
    else
      i = p:-1:1;
    end
    [alpha mu Xr] = varbvsbinupdate(X,sa,logodds,stats,alpha,mu,Xr,i);

    % Recalculate the posterior variance of the coefficients.
    s = sa./(sa*stats.xdx + 1);  

    % (2b) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    % -----------------------------------------------------
    % Compute the maximum a posteriori estimate of sa, if requested. Note
    % that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated. 
    if update_sa
      sa = (sa0*n0 + dot(alpha,s + mu.^2))/(n0 + sigma*sum(alpha));
      s  = sa./(sa*stats.xdx + 1);
    end
    
    % (2c) UPDATE ETA
    % ---------------
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    if update_eta
      eta   = update_eta(X,y,betavar(alpha,mu,s),Xr,stats.d);
      stats = update_stats(X,y,eta);
    end

    % (2d) COMPUTE VARIATIONAL LOWER BOUND
    % ------------------------------------
    % Compute variational lower bound to marginal log-likelihood.
    lnZ = intlogit(y,stats,alpha,mu,s,Xr,eta) ...
          + intgamma(logodds,alpha) ...
          + intklbeta(alpha,mu,s,sa);
    
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
     fprintf('%4d %+13.6e %0.1e %4d %0.2f\n',iter,lnZ,max(err),...
	      round(sum(alpha)),max(abs(alpha.*mu)));
    end
    if logw < logw0
      alpha = alpha0;
      mu    = mu0;
      eta   = eta0;
      logw  = logw0;
      break
    elseif max(err) < tol
      break
    end
  end

% ----------------------------------------------------------------------
% Calculates useful quantities for updating the variational approximation
% to the logistic regression factors.
function stats = update_stats (X, y, eta)

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute beta0 and yhat. See the journal paper for an explanation of
  % these two variables.
  beta0 = sum(y - 0.5)/sum(d);
  yhat  = y - 0.5 - beta0*d;

  % Calculate xy = X'*yhat as (yhat'*X)' and xd = X'*d as (d'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xd = double(d'*X)';

  % Compute the diagonal entries of X'*dhat*X. For a definition of dhat, see
  % the Bayesian Analysis journal paper.
  xdx = diagsq(X,d) - xd.^2/sum(d);

  % Return the result.
  stats = struct('d',d,'yhat',yhat,'xy',xy,'xd',xd,'xdx',xdx);

% ----------------------------------------------------------------------
% UPDATE_ETA(X,Y,V,XR,D) returns the M-step update for the parameters
% specifying the variational lower bound to the logistic regression factors.
%
% The inputs are as follows. Input X is an N x P matrix of observations
% about the variables (or features), where N is the number of samples, and P
% is the number of variables. Y is the vector of observations about the
% binary trait; it is a vector of length N. Vector Y and the columns of
% matrix X must not be centered. Input XR must be equal to XR = X*R, where R
% is the posterior mean of the coefficients. Note that under the
% fully-factorized variational approximation, R = ALPHA.*MU; see VARBVSBIN
% for details.
%
% Input V is the posterior variance of the coefficients. For this update to
% be valid, it is required that the posterior covariance of the coefficients
% is equal to DIAG(V). Input D must be D = SLOPE(ETA), where ETA is the
% current estimate of the free parameters; see function SLOPE for details.
function eta = update_eta (X, y, v, Xr, d)
  
  % Compute MU0, the posterior mean of the intercept in the logistic
  % regression under the variational approximation. Here, A is the variance
  % of the intercept given the other coefficients.
  a   = 1/sum(d);
  mu0 = a*(sum(y - 0.5) - d'*Xr);

  % Compute S0, the (marginal) posterior variance of the intercept in the
  % logistic regression. Here, I calculate XD = X'*D as (D'*X)' to avoid
  % storing the transpose of X, since X may be large.
  xd = double(d'*X)';
  s0 = a*(1 + a*v'*(xd.^2));
  
  % Calculate the covariance between the intercept and coefficients.
  c = -a*xd.*v;
  
  % This is the M-step update for the free parameters.
  eta = sqrt((mu0 + Xr).^2 + s0 + diagsqt(X,v) + 2*double(X*c));

% ----------------------------------------------------------------------
% INTLOGIT(Y,STATS,ALPHA,MU,S,XR,ETA) computes an integral that appears in
% the variational lower bound of the marginal log-likelihood for the
% logistic regression model. This integral is an approximation to the
% expectation of the logistic regression log-likelihood taken with respect
% to the variational approximation. 
%
% Input arguments Y and STATS specify the variational approximation to the
% likelihood. Y is the column vector of observations about the binary
% outcome. STATS is the STRUCT output from UPDATESTATS. Inputs ALPHA, MU and
% S are the parameters of the approximating distribution; see VARBVSBIN for
% details. ETA is the vector of free parameters specifying the variational
% lower bound to the logistic regression factors. Input XR must be equal to
% XR = X*R, where R is the posterior mean of the coefficients. Note that
% under the fully-factorized variational approximation, R = ALPHA.*MU; see
% VARBVSBIN for details. ETA and XR are column vectors with the same length
% as Y.
function I = intlogit (y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  xdx  = stats.xdx;
  d    = stats.d;
  D    = diag(sparse(d));

  % Get the variance of the intercept given the other coefficients.
  a = 1/sum(d);

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the variational approximation.
  I = sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 + log(a)/2 ...
      + a*sum(y - 0.5)^2/2 + yhat'*Xr - qnorm(Xr,D)^2/2 ...
      + a*(d'*Xr)^2/2 - xdx'*betavar(alpha,mu,s)/2;
