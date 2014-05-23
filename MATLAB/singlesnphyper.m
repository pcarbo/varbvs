% [ALPHA,MU,S] = SINGLESNPHYPER(X,Y,W,LOG10SIGMA,H,THETA0) computes
% posterior probabilities and expectations of the coefficients in Bayesian
% variable selection averaged over hyperparameter settings, assuming that
% all variables are independent of each other. Input W is the array of
% normalized importance weights. For an explanation of the other inputs,
% and the outputs, see function MULTISNPHYPER.
function [alpha, mu, s] = singlesnphyper (X, y, w, log10sigma, h, theta0)

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
  for i = 1:ns

    % Calculate the mean (mu) and variance of the coefficients (s) given
    % that the coefficients are included in the model, and the posterior
    % inclusion probabilities (alpha), ignoring correlations between
    % markers.
    [alpha(:,i) mu(:,i) s(:,i)] = indbvs(xy,d,sigma(i),sa(i),...
                                         log(10)*theta0(i));
  end

  % Compute the posterior expectations averaged over the settings of the
  % hyperparameters.
  alpha = alpha * w(:);
  mu    = mu    * w(:);
  s     = s     * w(:);
