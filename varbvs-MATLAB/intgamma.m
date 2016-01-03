% intgamma(logodds,alpha) computes an integral that appears in the
% variational lower bound of the marginal log-likelihood. This integral is
% the expectation on the prior inclusion probabilities taken with respect to
% the variational approximation. 
function y = intgamma (logodds, alpha)

  % This is the same as 
  %
  %    sum(alpha.*log(q) + (1-alpha).*log(1-q)).
  %  
  y = sum((alpha-1).*logodds + logsigmoid(logodds));
