% [ALPHA,MU,S] = INDBVS(XY,D,SIGMA,SA,LOGODDS) computes posterior
% probabilities and expectations of the coefficients in Bayesian variable
% selection assuming that all variables are independent of each other. See
% VARBVS and VARBVSUPDATE for details about the inputs and outputs.
function [alpha, mu, s] = indbvs (xy, d, sigma, sa, logodds)

  % Calculate the mean (mu) and variance of the coefficients (s) given that
  % the coefficients are included in the model, then calculate the posterior
  % inclusion probabilities (alpha), ignoring correlations between all
  % variables (or features).
  s     = sa*sigma./(sa*d + 1);
  mu    = s.*xy/sigma;
  SSR   = mu.^2./s;
  alpha = sigmoid(logodds + (log(s/(sa*sigma)) + SSR)/2);

