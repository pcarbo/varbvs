% [A,MU,XR] = MULTISNPBINUPDATE(X,SB,LOGODDS,STATS,A0,MU0,XR0,SNPS) runs
% a single iteration of the coordinate ascent algorithm for optimizing
% the variational lower bound of the variable selection model for a
% binary trait (such as case-control status). 
%
% X specifies the genotype data. It is an N x P matrix, where N is the
% number of samples (individuals), and P is the number of variables (genetic
% loci, or SNPs). X must be SINGLE. Unlike MULTISNPUPDATE, the matrix X
% should not be centered.
%
% Note that this routine is implemented with the assumption that the data X
% is single floating-point precision (type HELP SINGLE), as opposed to the
% MATLAB default of double precision. This is useful for large data sets,
% because single precision requires half of the number of bits as double
% floating-point precision. If X is provided in another numerical
% representation, an error is reported.
%
% For information on input STATS, see the help for function UPDATESTATS. For
% details on the outputs and the remaining inputs, see the help for function
% MULTISNPUPDATE. Note that a residual variance parameter SIGMA is not
% needed to model a binary trait.
%
% The computational complexity is O(N*M), where M = LENGTH(SNPS).

