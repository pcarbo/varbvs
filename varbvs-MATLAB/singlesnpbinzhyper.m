% [ALPHA,MU,S] = SINGLESNPBINZHYPER(X,Z,Y,H,THETA0,ETA) computes posterior
% probabilities and expectations of the coefficients for Bayesian variable
% selection in logistic regression, allowing for covariates, assuming that
% all variables are independent of each other. The posterior quantities are
% computed for each setting of the hyperparameters. For an explanation of
% the inputs and outputs, see function MULTISNPBINZHYPER.
function [alpha, mu, s] = singlesnpbinzhyper (X, Z, y, h, theta0, eta)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);

  % Get the settings for the prior variance of the additive effects.
  sx = sum(var1(X));
  sa = pve2sa(sx,h,theta0);

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
    fprintf('(%03d) h = %0.3f, theta0 = %+0.2f (sd = %0.3f)',...
	    i,h(i),theta0(i),sqrt(sa(i)));
    fprintf(repmat('\b',1,44));

    % Calculate the mean (mu) and variance of the coefficients (s) given
    % that the coefficients are included in the model, and the posterior
    % inclusion probabilities (alpha), ignoring correlations between
    % markers.
    [alpha(:,i) mu(:,i) s(:,i)] = ...
        indbvszbin(X,Z,y,sa(i),log(10)*theta0(i),eta(:,i));
  end
  fprintf('\n');
