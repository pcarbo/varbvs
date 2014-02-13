% *** DESCRIPTION OF FUNCTION GOES HERE ***
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
  