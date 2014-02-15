<% PVE(SX,SA,THETA0) returns the prior estimate of the proportion of variance
% in a given outcome (e.g. a quantitative trait) that is explained by a
% collection of explanatory variables (e.g. genetic markers). SA is the
% prior variance of the additive effects on the outcome, THETA0 is the
% (base 10) logarithm of the prior odds for inclusion. SX is the sum of the
% sample variances for all the explanatory variables. If inputs SA and
% THETA0 are both not scalars, they must be numeric arrays of the same
% dimension.
function h = pve (sx, sa, THETA0)

  % This is the expected value of the sample genetic variance divided by
  % the variance of the residual.
  c = sa .* sigmoid10(theta0) * sx;

  % This is the prior estimate of the proportion of variance in the outcome
  % that is explained by the variables. 
  h = c./(c + 1);
