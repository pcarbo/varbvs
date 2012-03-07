% INTKLBETA(A,MU,S,SB) computes the negative Kullback-Leibler divergence
% between the variational approximation and the prior on the regression
% coefficients. Input SB specifies the prior variance of the regression
% coefficients times the residual variance (for binary traits, set the
% residual variance to 1). 
%
% Inputs A, MU and S are each vectors of the same length specifying the
% variational approximation; the Kth regression coefficient is normal with
% probability A(K), and zero with probability 1 - A(K). MU(K) and S(K) are
% the mean and variance of the normal mixture component. 
function I = intklbeta (alpha, mu, s, sb)
  I = (sum(alpha) + alpha'*log(s/sb) - alpha'*(s + mu.^2)/sb)/2 ...
      - alpha'*log(alpha + eps) - (1 - alpha)'*log(1 - alpha + eps);
  