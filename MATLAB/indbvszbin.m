% DESCRIPTION OF FUNCTION GOES HERE.
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
  mu    = s.*xy/sigma;
  SSR   = mu.^2./s;
  alpha = sigmoid(logodds + (log(s/sa) + SSR)/2);
