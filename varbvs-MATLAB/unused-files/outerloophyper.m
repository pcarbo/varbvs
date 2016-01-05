% [LOGW,SIGMA,ALPHA,MU,S] = OUTERLOOPHYPER(X,Y,ALPHA,MU,SIGMA,H,THETA0,...
% VERBOSE) computes unnormalized log-importance weights for the
% hyperparameters. It is used by MULTISNPHYPER to implement the "outer loop"
% of the inference algorithm for analysis of a quantitative trait.
function [logw, sigma, alpha, mu, s] = ...
    outerloophyper (X, y, alpha, mu, sigma, h, theta0, verbose)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);

  % Get the settings for the prior variance of the additive effects.
  sx = sum(var1(X));
  sa = pve2sa(sx,h,theta0);

  % Repeat for each combination of the hyperparameters.
  for i = 1:ns
    fprintf('(%03d) h = %0.3f, theta = %+0.2f (sd = %0.3f) ',...
            i,h(i),theta0(i),sqrt(sa(i)));

    % Compute the unnormalized log-importance weight given values for the
    % hyperparameters, LOG10SIGMA, H and THETA0. Implicitly, the importance
    % weight includes these terms: the likelihood, the prior, and the
    % proposal. The proposal and prior cancel out from the expression for
    % the importance weight because both are assumed to be uniform for all
    % the hyperparameters.
    options = struct('alpha',alpha(:,i),'mu',mu(:,i),'update_sigma',true,...
                     'verbose',verbose);
    [logw(i) alpha(:,i) mu(:,i) s(:,i) sigma(i)] = ...
	varbvs(X,y,sigma(i),sa(i),log(10)*theta0(i),options);
  end
