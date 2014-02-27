% [X,Y] = CREATEBINDATA(MAF,BETA,MU,N) generates N samples of the genotypes
% and the binary outcome Y, according to SNP minor allele frequencies MAF,
% intercept MU, and log-odds ratios BETA. Inputs MAF and BETA are vectors of
% length P. The genotype data X is an N x P matrix, and the binary trait
% data Y is a column vector of length N. Note that, unlike function
% CREATEDATA, the columns of X is are not centred.
%
% Genotypes are generated from a binomial distribution with success rates
% given by the minor allele frequencies. Observations about the binary trait
% are generated according to a logistic regression model with log-odds of
% success (Y = 1) given by MU + X*BETA.
function [X, y] = createbindata (maf, beta, mu, n)

  % Get the number of SNPs.
  p = length(maf);
  
  % Simulate genotype data X assuming all the SNPs are uncorrelated (i.e. no
  % linkage disequilibrium between the SNPs), according to the specified
  % minor allele frequencies (MAFs).
  X = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));

  % For each subject, calculate the probability of having the disease (y = 1).
  p1 = sigmoid(mu + X*beta);

  % Simulate the binary trait (case-control status) as a coin toss with
  % success rates given by the logistic regression.
  y = double(rand(n,1) < p1);
  