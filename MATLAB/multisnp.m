% [LNZ,ALPHA,MU,S] = MULTISNP(X,Y,SIGMA,SB,LOGODDS,OPTIONS) finds the best
% variational approximation to the variable selection model for a
% quantitative trait.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (genetic loci,
% or SNPs). Y is the vector of quantitative trait data; it is a vector of
% length N. Crucially, this algorithm will only work correctly if Y and X
% are centered so that vector Y and each column of X has a mean of zero.
%
% Inputs SIGMA and SB are positive scalars. SIGMA specifies the variance of
% the residual, and SB is the prior variance of the additive effects.
% LOGODDS is the prior log-odds for each SNP. It is equal to LOGODDS =
% LOG(Q./(1-Q)), where Q is the prior probability that each SNP is included
% in the linear model of Y. LOGODDS may either be a scalar, in which case
% all the SNPs have the same prior, or it may be a vector of length P.
%
% Output LNZ is the variational estimate of the marginal log-likelihood
% given the hyperparameters SIGMA, SB and LOGODDS.
%
% Outputs ALPHA, MU and S are the parameters of the variational
% approximation; the Kth regression coefficient is normal with probability
% A(K), and zero with probability 1 - A(K). MU(K) and S(K) are the mean and
% variance of the normal mixture component.
%
% OPTIONS is an optional argument that overrides the default behaviour of
% the algorithm. Set OPTIONS.ALPHA and OPTIONS.MU to override the random
% initialization of the variational parameters ALPHA and MU. And set
% OPTIONS.VERBOSE = FALSE to turn off reporting progress of the algorithm.
function [lnZ, alpha, mu, s] = multisnp (X, y, sigma, sb, logodds, options)

  % Convergence criterion for coordinate ascent.
  tolerance = 1e-4;  
  
  % Get the number of participants in the study (n) and the number of SNPs
  % genotyped (p).
  [n p] = size(X);

  % TAKE CARE OF OPTIONS.
  if ~exist('options')
    options = [];
  end

  % Set initial estimates of sufficient statistics.
  if isfield(options,'alpha') & isfield(options,'mu')
    alpha = options.alpha;
    mu    = options.mu;
  else
    
    % Set the initial estimates of the sufficient statistics of the additive
    % effects so that exactly one SNP explains disease risk.
    alpha = rand(p,1);
    alpha = alpha / sum(alpha);
    mu    = randn(p,1);
  end

  % Determine whether to display the algorithm's progress.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end
  clear options
  
  % Make sure the genotypes are in single precision.
  if ~isa(X,'single')
    X = single(X);
  end

  % The prior log-odds should be a vector of length p.
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  else
    logodds = logodds(:);
  end

  % Compute a few useful quantities. Here, I calculate X'*Y as (Y'*X)' to
  % avoid storing the transpose of X, since X may be large. By keeping
  % track of the matrix-vector product X*r, we avoid having to recompute it
  % sometimes.
  xy = double(y'*X)';
  d  = diagsq(X);
  Xr = double(X*(alpha.*mu));
  
  % Calculate the variance of the additive effects.
  s = sb*sigma./(sb*d + 1);
  
  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  if verbose
    fprintf('iter           lnZ max chg #snp  beta\n');
  end
  while true
    iter = iter + 1;
    
    % Save the current variational parameters, and the variational lower
    % bound.
    alpha0  = alpha;
    mu0     = mu;
    lnZold  = lnZ;
    params0 = [ alpha; alpha .* mu ];

    % UPDATE VARIATIONAL APPROXIMATION.
    % Run a single forward or backward pass of the coordinate ascent
    % updates for the additive effects at each locus.
    if isodd(iter)
      snps = 1:p;
    else
      snps = p:-1:1;
    end
    [alpha mu Xr] = multisnpupdate(X,sigma,sb,logodds,xy,d,alpha,mu,Xr,snps);
    
    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute variational lower bound to the marginal log-likelihood.
    lnZ = intlinear(Xr,d,y,sigma,alpha,mu,s) ...
	  + intgamma(logodds,alpha) ...
	  + intklbeta(alpha,mu,s,sigma*sb);
    
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
    if lnZ < lnZold
      alpha = alpha0;
      mu    = mu0;
      lnZ   = lnZold;
      break
    elseif max(err) < tolerance
      break
    end
  end
