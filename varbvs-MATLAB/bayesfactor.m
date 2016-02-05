%--------------------------------------------------------------------------
% bayesfactor.m: Compute numerical estimate of Bayes factor.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Computes numerical estimate of
%
%           Pr(data | H1)
%      BF = ------------- ,
%           Pr(data | H0)
%
%    the probability of the data given the alternative hypothesis (H1) over
%    the probability of the data given the null hypothesis (H0). This is
%    also known as a Bayes factor (see Kass & Raftery, Journal of the
%    American Statistical Association, 1995). Here we assume that although
%    these probabilities cannot be computed analytically because they
%    involve intractable integrals, we can obtain reasonable estimates of
%    these probabilities with a simple numerical approximation over some
%    latent variable, Z. The inputs are the log-probabilities
%    Pr(data, Z | H0) and Pr(data, Z | H1) at different settings of Z.
%
% USAGE:
%    BF = varbvsbayesfactor(logw0, logw1)
%
% INPUT ARGUMENTS:
% logw0  log-probabilities under H0.
% logw1  log-probabilities under H1.
%
% OUTPUT ARGUMENTS: The estimated Bayes factor.
%
% LICENSE: GPL v3
%
% DATE: February 3, 2016
%
% AUTHORS:
%    Algorithm was designed by Peter Carbonetto and Matthew Stephens.
%    R, MATLAB and C code was written by Peter Carbonetto.
%    Depts. of Statistics and Human Genetics, University of Chicago,
%    Chicago, IL, USA, and AncestryDNA, San Francisco, CA, USA
%
% REFERENCES:
%    P. Carbonetto, M. Stephens (2012). Scalable variational inference
%    for Bayesian variable selection in regression, and its accuracy in 
%    genetic association studies. Bayesian Analysis 7: 73-108.
%
% SEE ALSO:
%    varbvs
%
function BF = varbvsbayesfactor (logw0, logw1)

  % Compute the marginal log-likelihood under the null hypothesis.
  c     = max(logw0(:));
  logz0 = c + log(mean(exp(logw0(:) - c)));

  % Compute the marginal log-likelihood under the alternative hypothesis.
  c     = max(logw1(:));
  logz1 = c + log(mean(exp(logw1(:) - c)));

  % Compute the numerical estimate of the Bayes factor.
  BF = exp(logz1 - logz0);
