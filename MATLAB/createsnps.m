% [MAF,BETA] = CREATESNPS(P,N) generates minor allele frequencies (MAF) and
% additive effects (BETA) for P genetic loci (specifically, these are single
% nucleotide polymorphisms, or SNPs for short). Additive effects are
% generated from the standard normal. Minor allele frequencies are uniform
% between 0.05 and 0.5. Input argument N specifies the number of causal
% SNPs (SNPs that have a nonzero additive effect on the trait).
function [maf, beta] = createsnps (p, n)
  
  % Generate additive effects for the SNPs, such that N of them have a
  % nonzero effect on the trait.
  I       = randperm(p);
  I       = I(1:n);
  beta    = zeros(p,1);
  beta(I) = randn(n,1);
  
  % Generate the minor allele frequencies. They are uniform on [0.05,0.5].
  maf = 0.05 + 0.45 * rand(1,p);
  