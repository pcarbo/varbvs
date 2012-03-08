% *** FIX THESE COMMENTS ***
% 
% [SIGMA,SB,Q,PIP,MU] = VSMCMC(X,Y,A,B,LOGP,M0,NS) simulates from the
% posterior distribution of the Bayesian variable selection model for a
% quantitative trait, using a Markov chain Monte Carlo method. This
% function is used in the simulation studies of the Bayesian Analysis paper.
% 
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (genetic loci,
% or SNPs). Y is the vector of quantitative trait data; it is a vector of
% length N. To account for an intercept, center input vector Y and input
% matrix X so that Y and each column of X has a mean of zero.
% 
% Input arguments A, B and LOGP specify the priors on the hyperparameters.
% Positive scalars A and B are the prior sample sizes for the Beta prior on
% the prior probability of inclusion, Q. LOGP is the logarithm of the prior
% density for variance parameter SB. It is a handle to a function that takes
% as input a single argument, and returns the log-density at that point. See
% FUNCTION_HANDLE for more information. We assume the standard
% noninformative prior on the residual variance parameter SIGMA.
% 
% M0 determines the maximum number of SNPs that can be included in the
% model. NS specifies the length of the Markov chain, and the number of
% Monte Carlo samples generated.
%
% Outputs SIGMA, SB and Q are simulated samples of the hyperparameters. They
% are all vectors of length NS. Output PIP is the vector of Monte Carlo
% estimates of the posterior inclusion probabilities for all the SNPs;
% PIP(K) is the posterior probability that the Kth SNP is included in the
% model of Y. And MU is the vector of mean additive effects for the SNPs.
function [SIGMA, SA, Q, PIP, mu] = vsmcmc (X, y, a, b, logp, m0, ns)

  % These two parameters specify the inverse gamma prior on the variance of
  % the residual (sigma). I use the uninformative prior, which is the
  % limit of the inverse gamma as the scale and shape parameters approach
  % zero.
  as = 0.01;
  bs = 0.01;
  
  % Get the number of samples (n) and the number of variables (p).
  [n p] = size(X);
  
  % Compute a couple useful quantities.
  yy = y'*y;
  xy = X'*y;

  % Initialize the Markov chain.
  sa    = 1;           % Prior variance of coefficients.
  gamma = zeros(p,1);  % Indicator variables.
  m     = 0;           % Number of included variables.
  g     = [];          % Set of included variables.
  S     = [];          % Inverse of X'*X.  

  % Initialize storage for the Markov chain.
  Q     = zeros(1,ns);  % Prior inclusion probability.
  PIP   = zeros(p,1);   % Posterior inclusion probabilities.
  mu    = zeros(p,1);   % Posterior means of coefficients.
  SIGMA = zeros(1,ns);  % Residual variance.
  SA    = zeros(1,ns);  % Prior variance of coefficients.
  
  % Repeat for each step in the Markov chain.
  for i = 1:ns
    fprintf('%6d (%03d)',i,m);
    fprintf(repmat('\b',1,12));
    
    % (1.) Sample sigma | gamma, sa.
    sigma = samplesigma(n,S,xy(g),yy,as,bs);

    % (2.) Sample sa | gamma, sigma.
    [sa S] = samplesa(X(:,g),xy(g),sigma,sa,logp);
    
    % (3.) Sample gamma | sigma, sa.
    [gamma g S] = samplegamma(X,S,xy,gamma,g,sigma,sa,a,b,m0);
    m = length(g);
    
    % Update the posterior inclusion probabilities and the posterior mean
    % coefficients. 
    if m > 0
      PIP   = PIP + gamma;
      mu(g) = mu(g) + S*xy(g);
    end
    
    % Store the samples of the hyperparameters.
    SIGMA(i) = sigma;
    SA(i)    = sa;
    Q(i)     = beta_rnd(1,a+m,b+p-m);
  end
  fprintf('\n');

  % Normalize the Monte Carlo estimates.
  PIP = PIP / ns;
  mu  = mu / ns;

% ------------------------------------------------------------------
% Compute the softmax function in a numerically stable manner.
function p = softmax (a)
  amax = max(a);
  p    = exp(a - amax);
  p    = p / sum(p);

% ------------------------------------------------------------------
% Compute the ratio of the sum of the exponentials in a numerically
% stable manner. The issue is that one of the numbers a or b could be
% very large. 
function y = ratiosumexp (a, b)
  amax = max(a);
  bmax = max(b);
  y    = exp(amax - bmax) * sum(exp(a - amax)) / sum(exp(b - bmax));
  
% ------------------------------------------------------------------
% Sample the residual variance (sigma) given the prior variance (sa) and
% the indicator variables (gamma). Note that sampling X from the inverse
% gamma distribution with shape A and scale B is equivalent to sampling
% 1/X from the gamma distribution with shape A and scale 1/B.
function sigma = samplesigma (n, S, xy, yy, as, bs)
  SSR   = dot(xy,S*xy);
  a     = (as + n)/2;         % The shape parameter.
  b     = 2/(bs + yy - SSR);  % The scale parameter.
  sigma = 1/(b*gamm_rnd(1,a));

% ------------------------------------------------------------------
% Sample the prior variance (sa) given the residual variance (sigma) and
% the indicator variables (gamma).
function [sa, S] = samplesa (X, xy, sigma, sa, logp)

  % Get the size of the data set (n) and the number of active SNPs (m).
  [n m] = size(X);

  % Draw a candidate from the random walk proposal.
  sanew = sa * exp(randn);
  
  % Compute the Metropolis-Hastings acceptance probability. Note that we
  % need to multiply this term by sanew/sa due to the fact that we are
  % not conducting the random walk over sa directly. This expression
  % appears to be fairly stable numerically.
  I      = eye(m);
  S      = inv(X'*X + I/sa);
  Snew   = inv(X'*X + I/sanew);
  SSR    = dot(xy,S*xy);
  SSRnew = dot(xy,Snew*xy);
  alpha  = exp(logp(sanew) - logp(sa) ...
	       + (m*log(sa/sanew) + logdet(Snew) - logdet(S) ...
		  + SSRnew/sigma - SSR/sigma)/2) * sanew/sa;

  % Take a Metropolis-Hastings step.
  if rand < min(1,alpha)
    sa = sanew;
    S  = Snew;
  end

% ------------------------------------------------------------------
% Construct a Markov chain on the indicator variables (gamma) using a
% random-walk Metropolis-Hastings step conditioned on a value of the
% prior variance (sa) and the global scale (sigma).
function [gamma, g, S] = samplegamma (X, S, xy, gamma, g, sigma, sa, a, b, m0)

  % Get the size of the data set (n), the number of SNPs (p), and
  % the number of selected SNPs.
  [n p] = size(X);
  m     = length(g);
  
  % Treat the limiting case when there are no active covariates.
  if m == 0
    
    % BIRTH MOVE (m = 0).
    % Compute the activation probabilities, then choose which
    % regression coefficient to activate, then compute the acceptance
    % probability. 
    [logr Snew] = logactprob0(X,xy,sigma,sa,1:p);
    c = max(logr);
    r = exp(logr - c);
    i = randtable(r);
    if c > 100
      alpha = 1;
    else
      alpha = 0.5 * a/(b+p-1) * exp(c) * sum(r);
    end
    if rand < min(1,alpha)
      
      % Activate the ith regression coefficient.
      gamma(i) = 1;
      g        = i;
      S        = Snew(i);
    end
  else

    % Figure out the probability of a birth move.
    if m == m0
      qbirth = 0;
    else
      qbirth = 0.5;
    end
    
    % BIRTH MOVE (m > 0).
    if rand < qbirth
      
      % Compute the activation probabilities.
      ng         = find(gamma == 0);
      [logr t z] = logactprob(X,xy,S,sigma,sa,g,ng);
      
      % Choose which regression coefficient to activate.
      k    = randtable(softmax(logr));
      i    = ng(k);
      t    = t(:,k);
      z    = z(k);
      gnew = [ g; i ];
      Snew = [ S + z*t*t'  -z*t
	         -z*t'       z  ];
      
      % Compute the acceptance probability.
      if m == m0 - 1
	qdeath = 1;
      else
	qdeath = 0.5;
      end
      [logrnew ans ans] = logdeactprob(xy(gnew),Snew,sigma,sa);
      alpha = 2*qdeath * (a+m)/(b+p-m-1) ...
	      * ratiosumexp(logr,logrnew - logrnew(end));
      
      if rand < min(1,alpha)
	
	% Activate the ith regression cofficient.
	gamma(i) = 1;
	g        = gnew;
	S        = Snew;
      end
    else
      
      % DEATH MOVE (m = 1).
      % Treat the limiting case when there is exactly one active
      % covariate. 
      if m == 1
	
	% Compute the acceptance probability.
	[logrnew ans] = logactprob0(X,xy,sigma,sa,1:p);
	c     = max(logrnew);
	rnew  = exp(logrnew - c);
	alpha = 2 * (b+p-1)/a * exp(-c) / sum(rnew);	
	if rand < min(1,alpha)
	  
	  % Deactivate the regression coefficient.
	  gamma(g) = 0; 
	  g        = [];
	  S        = [];
	end
	
      % DEATH MOVE (m > 1).
      else
	
	% Compute the deactivation probabilities.
	[logr t z] = logdeactprob(xy(g),S,sigma,sa);
	
	% Choose which regression coefficient to activate.
	k    = randtable(softmax(logr));
	nk   = [ 1:k-1 k+1:m ];
	i    = g(k);
	t    = t(nk,k);
	z    = z(k);
	gnew = g(nk);
	Snew = S(nk,nk) - z*t*t';	
	
	% Compute the acceptance probability.
	[logrnew ans ans] = logactprob(X,xy,Snew,sigma,sa,gnew,...
				       [ find(gamma == 0); i ]);
	alpha = 0.5/(1-qbirth) * (b+p-m-2)/(a+m-1) * ...
		ratiosumexp(logr,logrnew - logrnew(end));

	if rand < min(1,alpha)
	  
	  % Deactivate the ith regression coefficient.
	  gamma(i) = 0; 
	  g        = gnew;
	  S        = Snew;
	end
      end
    end
  end

% ------------------------------------------------------------------
% Compute the activation probabilities in the limiting case when there
% are no active covariates.
function [logr, S] = logactprob0 (X, xy, sigma, sa, is)
  S    = 1./(1/sa + diagsq(X(:,is)))';
  SSR  = (xy(is).^2)' .* S;
  logr = log(S/sa)/2 + SSR/(2*sigma);

% ------------------------------------------------------------------
% Compute the logarithm of the activation probabilities.
function [logr, t, z] = logactprob (X, xy, S, sigma, sa, g, is)
  a    = X(:,g)'*X(:,is);
  t    = S*a;
  z    = 1./(1/sa + diagsq(X(:,is))' - sum(a.*(S*a),1));
  logr = log(z/sa)/2 + z.*(xy(is)' - xy(g)'*t).^2/(2*sigma);

% ------------------------------------------------------------------
% Compute the logarithm of the deactivation probabilities for *all* 
% the active covariates.
function [logr, t, z] = logdeactprob (xy, S, sigma, sa)
  m     = length(xy);
  z     = diag(S)';
  t     = -(S - diag(diag(S)))./(repmat(z,m,1) + eps);
  XY    = repmat(xy,1,m);
  XY    = XY - diag(diag(XY));
  logr  = -(log(z/sa) + z.*(xy' - sum(XY.*t)).^2/sigma)/2;
