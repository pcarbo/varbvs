% *** FIX THESE COMMENTS ***
% 
% [LNZ,ALPHA,MU,S,ETA] = MULTISNPBIN(X,Y,SB,LOGODDS,OPTIONS) finds the best
% variational approximation to the variable selection model for a binary
% trait.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (genetic loci,
% or SNPs). Y is the vector of binary trait data, such as case-control
% status; it is a vector of length N. Unlike function MULTISNP, X and Y
% should not be centered. Instead, we will account for the intercept as
% we update the variational approximation.
%
% Input SB is a positive scalar specifying the prior variance of the
% additive effects. LOGODDS is the prior log-odds for each SNP. It is equal
% to LOGODDS = LOG(Q./(1-Q)), where Q is the prior probability that each SNP
% is included in the linear model of Y. LOGODDS may either be a scalar, in
% which case all the SNPs have the same prior, or it may be a vector of
% length P. Note that a residual variance parameter SIGMA is not needed to
% model a binary trait.
%
% Output LNZ is the variational estimate of the marginal log-likelihood
% given the hyperparameters SB and LOGODDS.
%
% Outputs ALPHA, MU and S are the parameters of the variational
% approximation; the Kth regression coefficient is normal with probability
% A(K), and zero with probability 1 - A(K). MU(K) and S(K) are the mean and
% variance of the normal mixture component. ETA is the set of parameters
% specifying the variational approximation to the nonlinear terms in the
% logistic regression. It is a vector of length N.
%
% OPTIONS is an optional argument that overrides the default behaviour of
% the algorithm. Set OPTIONS.ALPHA and OPTIONS.MU to override the random
% initialization of the variational parameters ALPHA and MU. Set OPTIONS.ETA
% to override the default initialization of the free parameters for the
% variational lower bound on the nonlinear terms in the logistic regression.
% Set OPTIONS.FIXEDETA = TRUE to keep the variational approximation for the
% logistic regression constant. And set OPTIONS.VERBOSE = FALSE to turn off
% reporting progress of the algorithm.
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
    error('Inputs SIGMA and SA must be scalars');
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
  if isfield(options,'alpha') & isfield(options,'mu')
    alpha = options.alpha(:);
    mu    = options.mu(:);
  else
    
    % The variational parameters are initialized randomly so that exactly
    % one variable explains the outcome.
    alpha = rand(p,1);
    alpha = alpha / sum(alpha);
    mu    = randn(p,1);
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
  if isfield(options,'fixedeta')
    fixedeta = options.fixedeta;
  else
    fixedeta = false;
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

  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  % *** FIX THIS ***
  if verbose
    fprintf('iter           lnZ max chg #snp   LOR\n');
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
      I = 1:p;
    else
      I = p:-1:1;
    end
    % *** FIX THIS ***
    [alpha mu Xr] = multisnpbinupdate(X,sa,logodds,stats,alpha,mu,Xr,snps);

    % Recalculate the posterior variance of the coefficients.
    s = sb./(sb*stats.d + 1);  

    % UPDATE ETA.
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    % *** FIX THIS ***
    if ~fixedeta
      eta   = updateeta(X,y,betavar(alpha,mu,s),Xr,stats.u);
      stats = updatestats(X,y,eta);
    end

    % UPDATE VARIATIONAL LOWER BOUND.
    % Compute variational lower bound to marginal log-likelihood.
    % *** FIX THIS ***
    lnZ = intlogit(y,stats,alpha,mu,s,Xr,eta) ...
	  + intgamma(logodds,alpha) ...
	  + intklbeta(alpha,mu,s,sa);
    
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [ alpha; alpha .* mu ];
    is     = find(abs(params) > 1e-6);
    err    = relerr(params(is),params0(is));
    if verbose
      fprintf('%4d %+13.6e %0.1e %4d %0.3f\n',iter,lnZ,max(err),...
	      round(sum(alpha)),max(abs(alpha.*mu)));
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
