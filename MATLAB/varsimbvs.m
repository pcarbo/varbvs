% [W,ALPHA,MU] = VARSIMBVS(X,Y,SIGMA,SA,LOG10Q,A,B,C) runs the full
% inference procedure for Bayesian variable selection in linear
% regression. This is a special implementation of the variational inference
% procedure used in the two simulation studies for the Bayesian Analysis
% paper. The main distinguishing feature of this procedure is the choice of
% priors for the hyperparameters of the variable selection model. In
% addition, we also avoid erratic behaviour in the approximation by first
% searching for a good initialization of the variational parameters. This
% inference procedure involves an inner loop and an outer loop. The inner
% loop consists of running a coordinate ascent algorithm to tighten the
% variational lower bound given a setting of the hyperparameters (this inner
% loop is implemented in function VARBVS). The outer loop computes
% importance weights for all combinations of the hyperparameters.
%
% Input X is an N x P matrix of observations about the variables (or
% features), where N is the number of samples, and P is the number of
% variables. Y is the vector of observations about the outcome; it is a
% vector of length N. To account for an intercept, Y and X must be centered
% beforehand so that Y and each column of X has a mean of zero.
%
% Note that this routine is implemented with the assumption that the data X
% is single floating-point precision (type HELP SINGLE), as opposed to the
% MATLAB default of double precision. This is useful for large data sets,
% because single precision requires half of the number of bits as double
% floating-point precision. If X is provided in another numerical
% representation, it is immediately converted to SINGLE.
%
% Inputs SIGMA, SA and LOG10Q specify the hyperparameter settings. These
% inputs must be arrays of the same size. For each combination of the
% hyperparameters, we compute an importance weight, and store the result in
% output W. SIGMA is the residual variance, SIGMA.*SA is the prior variance
% of the regression coefficients, and LOG10Q is the (base 10) logarithm of
% the prior inclusion probability.
%
% Inputs A, B and C are positive scalars. A and B are the prior sample sizes
% for the Beta prior on the prior inclusion probability. We assume a uniform
% prior on the proportion of variance explained, except that we replace the
% prior inclusion probability the proportion of variance explained by a
% constant, C. This is done purely for convenience, so that hyperparameter
% SA does not depend on the prior inclusion probability a priori, making it
% easier to implement the Markov chain Monte Carlo (MCMC) method. We assume
% the standard noninformative prior on the residual variance SIGMA.
%
% Outputs ALPHA and MU are variational estimates of the posterior inclusion
% probabilities and posterior mean of the coefficients (given that the
% variable is included in the model). These variational estimates are
% averaged over the settings of the hyperparameters.
function [w, alpha, mu] = varsimbvs (X, y, sigma, sa, log10q, a, b, c)
  
  % These two parameters specify the inverse gamma prior on the variance of
  % the residual (sigma). I use the uninformative prior, which is the
  % limit of the inverse gamma as the scale and shape parameters approach
  % zero.
  as = 0.01;
  bs = 0.01;
  
  % Get the number of samples (n), the number of variables (p), and the
  % number of combinations of the hyperparameters (ns).
  [n p] = size(X);
  ns    = numel(sigma);

  % Inputs A, B and C must be scalars.
  if ~(isscalar(a) & isscalar(b) & isscalar(c))
    error('Inputs A, B and C must be scalars');
  end
 
  % Inputs SIGMA, SA and LOG10Q must have the same number of elements.
  if numel(sa) ~= ns | numel(log10q) ~= ns
    error('Inputs SIGMA, SA and LOG10Q must be the same size');
  end
  
  % Get the sum of the sample variances.
  sx = sum(var1(X));
  
  % Get the settings for the prior inclusion probabilities.
  q = 10.^log10q;

  % Initialize storage for the marginal log-likelihoods (lnZ), the
  % log-importance weights (logw), variational estimates of the posterior
  % inclusion probabilities (alpha), and variational estimates of the
  % posterior mean coefficients (mu).
  lnZ   = zeros(size(sa));
  logw  = zeros(size(sa));
  alpha = zeros(p,ns);
  mu    = zeros(p,ns);
  
  % First get the best initialization for the variational parameters. Repeat
  % for each combination of the hyperparameters.
  fprintf('Finding best initialization for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  for i = 1:ns
    fprintf('(%03d) sigma = %4.1f, sa = %0.3f, q = %0.2e',...
	    i,sigma(i),sa(i),q(i));
    fprintf(repmat('\b',1,44));
  
    % Randomly initialize the variational parameters.
    alpha0  = rand(p,1);
    alpha0  = alpha0 / sum(alpha0);
    options = struct('alpha',alpha0,'mu',randn(p,1),'verbose',false);

    % Run the coordinate ascent algorithm.
    [lnZ(i) alpha(:,i) mu(:,i)] = ...
	varbvs(X,y,sigma(i),sa(i),logit(q(i)),options);
  end
  fprintf('\n');
  
  % Choose an initialization common to all the runs of the coordinate ascent
  % algorithm. This is chosen from the hyperparameters with the highest
  % marginal likelihood.
  [ans i] = max(lnZ(:));
  options = struct('alpha',alpha(:,i),'mu',mu(:,i),'verbose',false);
  
  % Repeat for each combination of the hyperparameters.
  fprintf('Computing importance weights for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  for i = 1:ns
    fprintf('(%03d) sigma = %4.1f, sa = %0.3f, q = %0.2e',...
	    i,sigma(i),sa(i),q(i));
    fprintf(repmat('\b',1,44));
  
    % Run the coordinate ascent algorithm.
    [lnZ(i) alpha(:,i) mu(:,i)] = ...
	varbvs(X,y,sigma(i),sa(i),logit(q(i)),options);

    % Compute the log-importance weight. Note that if X ~ p(x), then the
    % probability density of Y = log(X) is proportional to x*p(x). This is
    % useful for calculating the prior for the logarithm of the prior
    % inclusion probability, since we want to calculate the posterior
    % distribution for log10(q), not q.
    logw(i) = lnZ(i) ...                            % Marginal likelihood.
	      + loginvgamma(sigma(i),as/2,bs/2) ... % Prior on sigma.
	      + logpve(c*sx,sa(i)) ...              % Prior on sa.
	      + logbeta(q(i),a+1,b);                % Prior on log10q.
  end
  fprintf('\n');
  
  % Compute the normalized importance weights.
  w = normalizelogweights(logw);

  % Compute the weighted averages.
  alpha = alpha * w(:);
  mu    = mu    * w(:);

