% [ALPHA,MU,S] = SINGLESNPHYPER(X,Y,LOG10SIGMA,H,THETA0) computes posterior
% probabilities and expectations of the coefficients for Bayesian variable
% selection in linear regression, assuming that all variables are
% independent of each other. The posterior quantities are computed for each
% setting of the hyperparameters. For an explanation of the inputs and
% outputs, see function MULTISNPHYPER.
function [alpha, mu, s] = singlesnphyper (X, y, log10sigma, h, theta0)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);

  % Get the settings for the prior residual variance (sigma), and the prior
  % variance of the additive effects (sa).
  sigma = 10.^log10sigma;
  sx    = sum(var1(X));
  sa    = pve2sa(sx,h,theta0);

  % Compute a couple useful quantities. Here I calculate X'*Y as (Y'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(y'*X)';
  d  = diagsq(X);

  % Initialize storage for the posterior inclusion probabilities (alpha),
  % posterior means (mu), and posterior variances (s) of the regression
  % coefficients.
  alpha = zeros(p,ns);
  mu    = zeros(p,ns);
  s     = zeros(p,ns);

  % Repeat for each combination of the hyperparameters.
  fprintf('Computing posterior expectations for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  for i = 1:ns
    fprintf('(%03d) sigma = %0.2f, h = %0.3f, theta = %+0.2f (sd = %0.3f)',...
	    i,sigma(i),h(i),theta0(i),sqrt(sa(i)));
    fprintf(repmat('\b',1,57));

    % Calculate the mean (mu) and variance of the coefficients (s) given
    % that the coefficients are included in the model, and the posterior
    % inclusion probabilities (alpha), ignoring correlations between
    % markers.
    [alpha(:,i) mu(:,i) s(:,i)] = ...
        indbvs(xy,d,sigma(i),sa(i),log(10)*theta0(i));
  end
  fprintf('\n');
