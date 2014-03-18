% PVE2SA(SX,H,THETA0) returns the prior variance of the additive effects on
% the outcome given the prior estimate of the proportion of variance
% explained by the additive effects (H). THETA0 is the (base 10) logarithm
% of the prior odds for inclusion, and SX is the sum of the sample variances
% for all the explanatory variables.
function sa = pve2sa (sx, h, theta0)
  sa = h./(1-h)./(sigmoid10(theta0) * sx);
