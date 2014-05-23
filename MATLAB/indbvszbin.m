% [ALPHA,MU,S] = INDBVSZBIN(X,Z,Y,SA,LOGODDS,ETA) computes posterior
% probabilities and expectations of the coefficients for Bayesian variable
% selection in logistic regression, allowing for covariates, ignoring
% correlations between the variables. See VARBVSZBIN and VARBVSZBINUPDATE
% for details about the inputs and outputs.
function [alpha, mu, s] = indbvszbin (X, Z, y, sa, logodds, eta)

  % Get the statistics used to calculate the posterior probabilities and
  % expectations below.
  stats = updatestats2(X,Z,y,eta);
  xy    = stats.xy;
  xdx   = stats.xdx;

  % Calculate the mean (mu) and variance of the coefficients (s) given that
  % the coefficients are included in the model, then calculate the posterior
  % inclusion probabilities (alpha), ignoring correlations between all
  % variables (or features).
  s     = sa./(sa*xdx + 1);
  mu    = s.*xy;
  SSR   = mu.^2./s;
  alpha = sigmoid(logodds + (log(s/sa) + SSR)/2);
