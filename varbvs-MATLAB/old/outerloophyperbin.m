% [LOGW,ALPHA,MU,S,ETA] = OUTERLOOPHYPERBIN(X,Y,ALPHA,MU,ETA,H,THETA0)
% computes unnormalized log-importance weights for the hyperparameters. It
% is used by MULTISNPBINHYPER to implement the "outer loop" of the inference
% algorithm for analysis of a binary trait.
function [logw, alpha, mu, s, eta] = ...
      outerloophyperbin (X, y, alpha, mu, eta, h, theta0)

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
    fprintf('(%03d) h = %0.3f, theta0 = %+0.2f (sd = %0.3f)\n',...
	    i,h(i),theta0(i),sqrt(sa(i)));
  
    % Compute the unnormalized log-importance weight given values for the
    % hyperparameters, H and THETA0. Implicitly, the importance weight
    % includes these terms: the likelihood, the prior, and the proposal. The
    % proposal and prior cancel out from the expression for the importance
    % weight because both are assumed to be uniform for all the
    % hyperparameters.
    options = struct('alpha',alpha(:,i),'mu',mu(:,i),'eta',eta(:,i));
    [logw(i) alpha(:,i) mu(:,i) s(:,i) eta(:,i)] = ...
	varbvsbin(X,y,sa(i),log(10)*theta0(i),options);
    fprintf('\n');
  end
