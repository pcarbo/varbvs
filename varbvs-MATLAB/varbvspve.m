% TO DO: Update these comments.
%
% NOTE: Assumes intercept.
%
% [HPE,DPE] = SAMPLEPVE(X,K,Y,ALPHA,MU,S,W,LOG10SIGMA,LOG10ODDS,H,D,NS)
% samples posterior estimates of (1) the proportion of variance explained
% (output HPE), and (2) the proportion of additive genetic variance due to
% QTL effects (output DPE), under the hybrid model. Outputs HPE and DPE
% contain NS samples of these posterior estimates.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (genetic loci,
% or SNPs). Y is the vector of quantitative trait data. It is a vector of
% length N. Crucially, this algorithm will only work correctly if Y and X
% are centered so that vector Y and each column of X has a mean of zero. 
% Input K = X*X'/P is the kinship matrix.
%
% The combinations of hyperparameter settings are given by inputs
% LOG10SIGMA, LOG10ODDS, H and D. W is the array of importance weights for
% each combination of the hyperparameters. These inputs must be arrays with
% the same number of elements. For further details on these inputs, see
% function HYBRID. Inputs ALPHA, MU and S are variational estimates of the
% posterior probabilities and posterior moments for each combination of the
% hyperparameters. Each of these inputs are a P x NP matrix, where NP is the
% number of combinations of hyperparameters.
function pve = varbvspve (X, fit, nr)

  % Take care of the option inputs.
  if nargin < 3
    nr = 1000;
  end

  % Get the number of variables.
  p = size(X,2);
    
  % Initialize storage for posterior estimates of the proportion of variance
  % explained.
  pve = zeros(nr,1);

  % Compute the normalized (approximate) importance weights.
  w = normalizelogweights(fit.logw);

  % For each sample, compute the proportion of variance explained.
  for i = 1:nr

    % Draw a hyperparameter setting from the posterior distribution.
    j = randtable(w);
    
    % Sample the region coefficients.
    b = fit.mu(:,j) + sqrt(fit.s(:,j)) .* randn(p,1);
    b = b .* (rand(p,1) < fit.alpha(:,j));

    % Compute the proportion of variance explained.
    sz     = var(X*b,1);
    pve(i) = sz/(sz + fit.sigma(j));
  end  
