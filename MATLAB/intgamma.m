% INTGAMMA(LOGODDS,A) computes the expectation of the prior on the indicator
% variables, where the expectation is taken with respect to the
% fully-factorized variational approximation. 
%
% Input LOGODDS is either a scalar, or a vector, specifying the prior
% log-odds. LOGODDS is equal to LOG(Q./(1-Q)), where Q is the prior
% probability that each SNP is included in the linear model of Y. 
%
% Input A specifies the mixture weights for the variational approximation;
% the Kth regression coefficient is normal with probability A(K), and zero
% with probability 1 - A(K).
function I = intgamma (logodds, alpha)

  % This is the same as 
  %
  %    SUM(ALPHA.*LOG(Q) + (1-ALPHA).*LOG(1-Q)).
  %  
  I = sum((alpha-1) .* logodds + logsigmoid(logodds));
