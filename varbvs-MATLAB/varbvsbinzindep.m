% TO DO: Add description of function (see varbvsbin.m and varbvsbinz.m).
function [alpha, mu, s] = varbvsbinzindep (X, Z, y, eta, sa, logodds)

  % Get the statistics used to calculate the posterior probabilities and
  % expectations below.
  stats = updatestats2(X,Z,y,eta);
  xy    = stats.xy;
  xdx   = stats.xdx;

  % Calculate the mean (mu) and variance (s) of the coefficients given that
  % they are included in the model, then calculate the posterior inclusion
  % probabilities (alpha), ignoring correlations between variables.

  s     = sa./(sa*xdx + 1);
  mu    = s.*xy;
  alpha = sigmoid(logodds + (log(s/sa) + mu.^2./s)/2);
