% [X,Y] = CREATEDATA(MAF,BETA,SIGMA,N) generates N samples of the genotypes
% and the quantitative trait (the continuous outcome Y), according to SNP
% minor allele frequencies MAF, and additive effects BETA. Inputs MAF and
% BETA are vectors of length P. The genotype data X is an N x P matrix, and
% the quantitative trait data Y is a column vector of length N. Both X and Y
% are centered so that Y and each column of X has a mean of zero.
%
% Genotypes are generated from a binomial distribution with success rates
% given by the minor allele frequencies. Observations about the quantitative
% trait are generated according to Y = X*BETA + E, where the residual E is
% normal with mean zero and covariance sigma*I.
function [X, y] = createdata (maf, beta, sigma, n)

  % Get the number of SNPs.
  p = length(maf);
  
  % Simulate genotype data X from an idealized population, according to the
  % specified minor allele frequencies.
  X = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));

  % Center the columns of X.
  X = X - repmat(mean(X),n,1);

  % Generate the quantitative trait measurements.
  y = X*beta + sqrt(sigma) * randn(n,1);

  % Take into account an intercept by centering the outcomes Y to have
  % mean zero.
  y = y - mean(y);
