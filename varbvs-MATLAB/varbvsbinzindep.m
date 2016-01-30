% This function computes the mean (mu) and variance (s) of the coefficients
% given that they are included in the logistic regression model, then it
% computes the posterior inclusion probabilities (alpha), ignoring
% correlations between variables. This function is used in varbvsindep.m.
function [alpha, mu, s] = varbvsbinzindep (X, Z, y, eta, sa, logodds)
  stats = update_varbvsbinz_stats(X,Z,y,eta);
  s     = sa./(sa*stats.xdx + 1);
  mu    = s.*stats.xy;
  alpha = sigmoid(logodds + (log(s/sa) + mu.^2./s)/2);
