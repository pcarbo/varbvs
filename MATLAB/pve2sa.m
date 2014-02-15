% PVE2SA(SX,H,LOG10ODDS) returns the prior variance of the additive effects
% on the outcome given the prior estiamte of the proportion of variance
% explained by the additive effects (H). LOG10ODDS is the (base 10)
% logarithm of the prior odds for inclusion, and SX is the sum of the sample
% variances for all the explanatory variables. If inputs SA and LOG10ODDS
% are both not scalars, they must be numeric arrays of the same dimension.
function sa = pve2sa (sx, h, log10odds)
  sa = h./(1-h)./(sigmoid10(log10odds) * sx);
