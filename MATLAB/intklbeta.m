% INTKLBETA(ALPHA,MU,S,SA) computes an integral that appears in the
% variational lower bound of the marginal log-likelihood. This integral is
% the negative Kullback-Leibler divergence between the approximating
% distribution and the prior of the coefficients. Input SA specifies the
% prior variance of the coefficients. (This SA is not the same as the SA
% used as an input to VARBVS.) See VARBVS for details on the inputs to this
% function.
function I = intklbeta (alpha, mu, s, sa)
  I = (sum(alpha) + alpha'*log(s/sa) - alpha'*(s + mu.^2)/sa)/2 ...
      - alpha'*log(alpha + eps) - (1 - alpha)'*log(1 - alpha + eps);
  