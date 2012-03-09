% [ALPHA,MU,XR] = VARBVSUPDATE(X,SIGMA,SA,LOGODDS,XY,D,ALPHA0,MU0,XR0,I)
% runs a single iteration of the coordinate ascent updates to maximize the
% variational lower bound for Bayesian variable selection in linear
% regression. It adjusts the fully-factorized variational approximation to
% the posterior distribution of the coefficients in a linear regression
% model of a continuous outcome (quantitiative trait), with spike and slab
% priors on the coefficients.
%
% All inputs to this function are required. Input X is an N x P matrix of
% observations about the variables (or features), where N is the number of
% samples, and P is the number of variables. X must be a SINGLE matrix.
% Input XY = X'*Y, where Y is the vector of observations about the outcome.
% Crucially, to account for an intercept, Y and X must be centered
% beforehand so that Y and each column of X has a mean of zero.
%
% Inputs SIGMA, SA and LOGODDS specify the hyperparameters. SIGMA and SA are
% scalars. SIGMA specifies the variance of the residual, and SA*SIGMA is the
% prior variance of the additive effects. LOGODDS is the prior log-odds of
% inclusion for each variable. It is equal to LOGODDS = LOG(Q./(1-Q)), where
% Q is the prior probability that each variable is included in the linear
% model of Y. LOGODDS is a vector of length P.
%
% Inputs ALPHA0, MU0 are the current parameters of the variational
% approximation; under the variational approximation, the ith regression
% coefficient is normal with probability ALPHA0(i), and MU0(i) is the mean
% of the coefficient given that it is included in the model. Inputs XR0 and
% D must be XR0 = X*(ALPHA0.*MU0) and D = DIAG(X'*X).
%
% Input I specifies the order in which the coordinates are updated. It may
% be a vector of any length. Each entry of I must be an integer between 1
% and P.
%
% There are three outputs. Output vectors ALPHA and MU are the updated
% variational parameters, and XR = X*(ALPHA.*MU). The computational
% complexity of VARBVSUPDATE is O(N*LENGTH(I)).

