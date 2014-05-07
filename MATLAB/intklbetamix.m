% INTKLBETAMIX(ALPHA,MU1,MU2,S1,S2,SA1,SA2) computes an integral that
% appears in the variational lower bound of the marginal log-likelihood for
% the mixture model. This integral is the negative Kullback-Leibler
% divergence between the approximating distribution and the prior of the
% coefficients. Inputs SA1 and SA2 specify the prior variance of the
% coefficients in the two normal mixture components. (These are not the
% same as SA1 and SA2 used as an input to VARBVSIX.)
function I = intklbetamix (alpha, mu1, mu2, s1, s2, sa1, sa2)
  a = alpha;
  I = (sum(a) + a'*log(s1/sa1) - a'*(s1 + mu1.^2)/sa1)/2 ...
      + (sum(1-a) + (1-a)'*log(s2/sa2) - (1-a)'*(s2 + mu2.^2)/sa2)/2...
      - a'*log(a + eps) - (1-a)'*log(1-a + eps);
  