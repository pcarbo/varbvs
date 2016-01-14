% int_klbeta(alpha,mu,s,sa) computes an integral that appears in the
% variational lower bound of the marginal log-likelihood. This integral is
% the negative K-L divergence between the approximating distribution and the
% prior of the coefficients. Note that this sa is not the same as the sa
% used as an input to varbvsnorm.
function I = int_klbeta (alpha, mu, s, sa)
  I = (sum(alpha) + alpha'*log(s/sa) - alpha'*(s + mu.^2)/sa)/2 ...
      - alpha'*log(alpha + eps) - (1 - alpha)'*log(1 - alpha + eps);
  