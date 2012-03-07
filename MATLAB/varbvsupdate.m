% [A,MU,XR] = MULTISNPUPDATE(X,SIGMA,SB,LOGODDS,XY,D,A0,MU0,XR0,SNPS) runs a
% single iteration of the coordinate ascent algorithm for optimizing the
% variational lower bound of the variable selection model for a quantitative
% trait.
%
% X specifies the genotype data. It is an N x P matrix, where N is the
% number of samples (individuals), and P is the number of variables (genetic
% loci, or SNPs). X must be SINGLE. XY is equal to X'*Y, where Y is the
% vector of quantitative trait data. D is equal to DIAG(X'*X). It is
% important to note that this function will only work correctly if Y and X
% are centered so that Y and each column of X has a mean of zero.
%
% Inputs SIGMA and SB are positive scalars. SIGMA specifies the variance of
% the residual, and SB is the prior variance of the additive
% effects. LOGODDS is the prior log-odds for each SNP. It is equal to
% LOGODDS = LOG(Q./(1-Q)), where Q is the prior probability that each SNP is
% included in the linear model of Y. LOGODDS is a vector of length P.
%
% Input SNPS specifies the order in which the coordinates are updated. It
% may be a vector of any length. Each entry of SNPs should be an integer
% between 1 and P.
%
% Inputs A0 and MU0 are vectors of length P specifying the initial
% variational approximation. Under the fully-factorized variational
% approximation, the Kth regression coefficient (additive effect) is normal
% with probability A0(K), and zero with probability 1 - A0(K). The mean of
% the normal is given by the Kth entry of MU0. Output vectors A and MU are
% the updated parameters of this variational approximation. Input XR0 must
% be equal to X*R0, such that R0 = A0.*MU0. Output XR is equal to X*R, where
% R = A.*MU.
%
% The computational complexity is O(N*M), where M = LENGTH(SNPS).
