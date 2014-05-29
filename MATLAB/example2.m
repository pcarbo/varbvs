% In this small example, I simulate a data set from an idealized genetic
% association study in which the genetic markers (SNPs) are uncorrelated. I
% assess the performance of the variational approximation for mapping
% associations between the genetic markers and a quantitative trait in this
% setting.
%
% Note that I have tested this script in MATLAB 8.1, and it may not work
% correctly in earlier versions of MATLAB.
clear

% SCRIPT PARAMETERS.
p  = 1e3;   % Number of variables (SNPs).
n  = 2e3;   % Number of samples (subjects).
na = 20;    % Number of SNPs that affect the outcome.
sb = 0.2;   % Standard deviation of nonzero coefficients.
se = 9;     % Variance of residual.

% Candidate values for the prior proportion of variance explained (h) and
% prior log-odds of inclusion (theta0).
theta0 = (-2.5:0.25:-1.5)';
h      = (0.05:0.05:0.45)';

% Set the random number generator seed.
seed = 1;
rng(seed);

% GENERATE DATA.
% Note that X and y are centered.
fprintf('Creating data set.\n');
[maf beta] = createsnps(p,na);
beta       = sqrt(se)*sb*beta;
[X y]      = createdata(maf,beta,se,n);

% Calculate the proportion of variance explained. Here, SZ is the sample
% genetic variance.
sz = var(X*beta,1);
fprintf('Proportion of variance explained is %0.3f.\n',sz/(sz + se));

% COMPUTE VARIATIONAL ESTIMATES.
fprintf('Computing variational estimates.\n');
[H THETA0] = ndgrid(h,theta0);
[logw sigma alpha mu s] = multisnphyper(X,y,H,THETA0);

% Compute the normalized importance weights.
w = normalizelogweights(logw);

% Show the posterior mean of each of the hyperparameters.
fprintf('Posterior mean of hyperparameters:\n');
fprintf('param   mean\n');
fprintf('sigma  %5.2f\n',dot(sigma(:),w(:)));
fprintf('theta0 %+0.2f\n',dot(THETA0(:),w(:)));
fprintf('h      %5.2f\n',dot(H(:),w(:)));
fprintf('\n');

% Calculate the posterior inclusion probabilities (PIPs), and show the PIPs
% for the variables that are most likely to be included in the model,
% with PIP > 0.1.
PIP     = alpha * w(:);
[ans I] = sort(-PIP);
I       = I(1:sum(PIP > 0.1));
fprintf('Selected variables:\n');
fprintf(' PIP    mu  beta\n');
fprintf('%0.2f %+0.2f %+0.2f\n',[PIP(I) mu(I) beta(I)]');
