%--------------------------------------------------------------------------
% varbvs.m: Fit variable selection model using variational approximation.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%   Compute fully-factorized variational approximation for Bayesian variable
%   selection in linear (family = 'gaussian') or logistic regression (family
%   = 'binomial'). More precisely, find the "best" fully-factorized
%   approximation to the posterior distribution of the coefficients, with
%   spike-and-slab priors on the coefficients. By "best", we mean the
%   approximating distribution that locally minimizes the Kullback-Leibler
%   divergence between the approximating distribution and the exact
%   posterior.
%
% USAGE:
%    fit = varbvs(X, Z, y)
%    fit = varbvs(X, Z, y, labels, family, options)
%    (Use empty matrix [] to apply the default value)
%  
% INPUT ARGUMENTS:
% X        n x p input matrix, where n is the number of samples,
%          and p is the number of variables. X cannot be sparse.
%
% Z        n x m covariate data matrix, where m is the number of
%          covariates. Do not supply an intercept as a covariate (i.e., a
%          column of ones), because an intercept is automatically included
%          in the regression model. For no covariates, set Z to the empty
%          matrix []. The covariates are assigned an improper, uniform
%          prior. Although improper priors are generally not advisable
%          because they can result in improper posteriors and Bayes factors,
%          this choice allows us to easily integrate out these covariates.
%
% y        Vector of length n containing observations of binary
%          (family = 'binomial') or continuous (family = 'gaussian')
%          outcome. For binary outcomes, all entries of y must be 
%          0 or 1.
%
% labels   Cell array with p entries containing variable labels. By
%          default, it is set to the empty matrix [].
%              
% family   Outcome type, either 'gaussian' or 'binomial'. Default
%          is 'gaussian'.
%
% options  A structure (type 'help struct') containing additional
%          model parameters and optimization settings. More details about
%          these options are given below. Fields, with their default
%          settings given, include:
%            
%          options.tol = 1e-4 (convergence tolerance for inner loop)
%          options.maxiter = 1e4 (maximum number of inner loop iterations)
%          options.verbose = true (print progress of algorithm to console)
%          options.sa = 1 (prior variance parameter settings)
%          options.logodds = linspace(-log10(p),-1,20) (log-odds settings)
%          options.update_sa (fit model parameter sa to data)
%          options.sa0 = 0 (scale parameter for prior on sa)
%          options.n0 = 0 (degrees of freedom for prior on sa)
%          options.alpha (initial estimate of variational parameter alpha)
%          options.mu (initial estimate of variational parameter mu)
%          options.initialize_params (see below)
%          options.nr = 100 (samples of PVE to draw from posterior)
%
%          For family = 'gaussian' only:
%          options.sigma = var(y) (residual variance parameter settings)
%          options.update_sigma (fit model parameter sigma to data)
%
%          For family = 'binomial' only:
%          options.eta (initial estimate of variational parameter eta)
%          options.optimize_eta (optimize parameter eta)
%
% OUTPUT ARGUMENTS:
% fit                 A structure (type 'help struct').
% fit.family          Either 'gaussian' or 'binomial'.
% fit.n               Number of training samples.
% fit.labels          Variable names.
% fit.sigma           settings for sigma (family = 'gaussian' only).
% fit.sa              settings for prior variance parameter, sa.
% fit.logodds         prior log-odds settings.
% fit.prior_same      true if prior is identical for all variables.
% fit.sa0             scale parameter for prior on sa.
% fit.n0              degrees of freedom for prior on sa.
% fit.update_sigma    whether sigma was fit to data (family = 'gaussian' only).
% fit.update_sa       whether hyperparameter sa was fit to data.
% fit.logw            approximate marginal log-likelihood for each
%                     setting of hyperparameters.
% fit.alpha           variational estimates of posterior inclusion probs.
% fit.mu              variational estimates of posterior mean coefficients.
% fit.s               variational estimates of posterior variances.
% fit.mu_cov          posterior estimates of coefficients for covariates.
% fit.eta             variational parameters for family = 'binomial' only.
% fit.optimize_eta    whether eta was fit to data (family = 'binomial' only).
% fit.pve             estimated PVE per variable (family = 'gaussian' only).
% fit.model_pve       posterior samples of proportion of variance in Y
%                     explained by variable selection model (only for
%                     family = 'gaussian').
%
% REGRESSION MODELS:
%    Two types of outcomes Y are modeled: (1) a continuous outcome, which
%    is also referred to as a "quantitative trait" in the genetics
%    literature; or (2) a binary outcome with possible values 0 and 1. Y is
%    modeled as a continuous outcome by setting family = 'gaussian'. In this
%    case, Y is i.i.d. normal with mean u0 + Z*u + X*b and variance sigma,
%    in which u and b are vectors of regresion coefficients, and u0 is the
%    intercept. In the second case, we use logistic regression to model Y,
%    in which the probability that Y = 1 is given by
%
%        Pr(Y = 1) = sigmoid(u0 + Z*u + X*b).
%
%    See 'help sigmoid' for a description of the sigmoid function. Note
%    that the regression always includes an intercept term (u0).
%
%    For both regression models, the fitting procedure consists of an inner
%    loop and an outer loop. The outer loop iterates over each of the
%    hyperparameter settings. The hyperparameters sa, sigma and logodds are
%    specified by three arrays with the same number of elements, in which
%    options.sa(i), options.sigma(i) and options.logodds(i) specify the ith
%    hyperparameter setting. (The exception to this is when logodds is a
%    matrix.) Note that sigma is only used for the linear regression model,
%    and will generate an error if family = 'binomial'.
%
% HYPERPARAMETERS:
%    Hyperparameter sa is the prior variance of regression coefficients for
%    variables that are included in the model. This prior variance is always
%    scaled by sigma (for logistic regression, we take sigma = 1).  Scaling
%    the variance of the coefficients in this way is necessary to ensure
%    that this prior is invariant to measurement scale (e.g., switching from
%    grams to kilograms). Hyperparameter logodds is the prior log-odds that
%    a variable is included in the regression model; it is defined as
%    logodds = log10(q/(1-q)), where q is the prior probability that a
%    variable is included in the regression model. Note that we use the
%    base-10 logarithm instead of the natural logarithm because it is
%    usually more natural to specify prior log-odds settings in this way.
%
%    The prior log-odds may also be specified separately for each variable,
%    which is useful is there is prior information about which variables are
%    most relevant to the outcome Y. This is accomplished by setting
%    options.logodds to a p x ns matrix, where p is the number of
%    variables, and ns is the number of hyperparameter settings. In this
%    case, fit.prior_same = false.
% 
% CO-ORDINATE ASCENT OPTIMIZATION:
%    Given a setting of the hyperparameters, options.sa(i), options.sigma(i)
%    and options.logodds(:,i), the inner loop cycles through coordinate
%    ascent updates to tighten the lower bound on the marginal likelihood,
%
%        Pr(Y | X, sigma, sa, logodds).
%
%    The inner loop coordinate ascent updates terminate when either (1) the
%    maximum number of inner loop iterations is reached, as specified by
%    options.maxiter, or (2) the maximum difference between the estimated
%    posterior inclusion probabilities (see below) is less than options.tol.
%
%    To provide a more accurate variational approximation of the posterior
%    distribution, by default the fitting procedure has two stages. In the
%    first stage, the entire fitting procedure is run to completion, and the
%    variational parameters (alpha, mu, s, eta) corresponding to the maximum
%    lower bound are then used to initialize the coordinate ascent updates
%    in a second stage. Although this has the effect of doubling the
%    computation time (in the worst case), the final posterior estimates
%    tend to be more accurate with this two-stage fitting procedure. The
%    initial stage will be automatically skipped if initial estimates of the
%    variational parameters are provided in options.alpha, options.mu and/or
%    options.s, unless options.initialize_params = true.
%
%    It is possible to optimize hyperparameter sa (and sigma for family =
%    'gaussian') as part of the inner loop fitting procedure. Parameters sa
%    and sigma will automatically be fitted to the data, separately for each
%    hyperparameter setting, when options.sa and options.sigma are not
%    specified. Alternatively, this can be achieved by setting
%    options.update_sa = true and options.update_sigma = true, in which case
%    options.sa and options.sigma are treated as initial estimates of these
%    parameters if they are provided. These parameters are fitted by
%    computing approximate maximum-likelihood (ML) estimates. Optionally, an
%    approximate maximum a posteriori (MAP) estimate of sa is computed by
%    setting options.sa0 and options.n0 to positive scalars; these two
%    numbers specify the scale parameter and number of degrees of freedom
%    for a scaled inverse chi-square prior on sa. Note it is not possible to
%    fit the logodds parameter; if options.logodds is not provided, then it
%    is set to the default value when options.sa and options.sigma are
%    scalars, and otherwise an error is generated.
%
% VARIATIONAL APPROXIMATION:
%    Outputs fit.alpha, fit.mu and fit.s specify the approximate posterior
%    distribution of the regression coefficients. Each of these outputs is a
%    p x ns matrix. For the ith hyperparameter setting, fit.alpha(:,i) is
%    the variational estimate of the posterior inclusion probability (PIP)
%    for each variable; fit.mu(:,i) is the variational estimate of the
%    posterior mean coefficient given that it is included in the model; and
%    fit.s(:,i) is the estimated posterior variance of the coefficient given
%    that it is included in the model. These are also the quantities that
%    are optimized as part of the inner loop coordinate ascent updates. From
%    these quantities, we also provide fit.pve for family = 'gaussian'. This
%    output provides, for each hyperparameter setting, the mean estimate of
%    the proportion of variance in Y explained by each of the variables
%    conditioned on being included in the model.
%
%    An additional variational parameter, denoted by 'eta', is needed for
%    fast computation with the logistic regression model (family =
%    'binomial'). The fitted value of eta is returned as an n x ns matrix
%    fit.eta. If a good estimate of eta is already available (e.g., in a
%    previous call to varbvs on the same data), provide this estimate in
%    options.eta, in which case eta is not fitted to the data during the
%    inner loop coordinate ascent updates (to override this behaviour, set
%    options.optimize_eta = true, in which case options.eta is treated as an
%    initial estimate).
%
%    The variational estimates should be interpreted carefully, especially
%    when variables are strongly correlated. For example, consider the
%    simple scenario in which 2 candidate variables are closely correlated,
%    and at least one of them explains the outcome with probability close to
%    1. Under the correct posterior distribution, we would expect that each
%    variable is included with probability ~0.5. However, the variational
%    approximation, due to the conditional independence assumption, will
%    typically get this wrong, and concentrate most of the posterior weight
%    on one variable (the actual variable that is chosen will depend on the
%    starting conditions of the optimization). Although the individual PIPs
%    are incorrect, a statistic summarizing the variable selection for both
%    correlated variables (e.g., the total number of variables included in
%    the model) should be reasonably accurate.
%
%    More generally, if variables can be reasonably grouped together based
%    on their correlations, we recommend interpreting the variable selection
%    results at a group level. For example, in genome-wide association
%    studies (see the vignettes) ,a SNP with a high PIP indicates that this
%    SNP is probably associated with the trait, and one or more nearby SNPs
%    within a chromosomal region, or ``locus,'' may be associated as
%    well. Therefore, we interpreted the GWAS variable selection results at
%    the level of loci, rather than at the level of individual SNPs.
%
%    Also note that special care is required for interpreting the results of
%    the variational approximation with the logistic regression model. In
%    particular, interpretation of the individual estimates of the
%    regression coefficients (e.g., the posterior mean estimates fit.mu) is
%    not straightforward due to the additional approximation introduced on
%    the individual nonlinear factors in the likelihood. As a general
%    guideline, only the relative magnitudes of the coefficients are
%    meaningful.
%
% AVERAGING OVER HYPERPARAMETER SETTINGS:
%    Output fit.logw is an array with ns elements, in which fit.logw(i) is
%    the variational lower bound on the marginal log-likelihood for the ith
%    setting of the hyperparameters. In many settings, it is good practice
%    to account for uncertainty in the hyperparameters when reporting final
%    posterior quantities. Provided that (1) the hyperparameter settings
%    options.sigma, options.sa and options.logodds adequately represent the
%    space of possible hyperparameter settings with high posterior mass, (2)
%    the hyperparameter settings are drawn from the same distribution as the
%    prior, and (3) the fully-factorized variational approximation closely
%    approximates the true posterior distribution, then final posterior
%    quantities can be calculated by using fit.logw as (unnormalized)
%    log-marginal probabilities. Even when conditions (1), (2) and/or (3)
%    are not satisfied, this can approach can still often yield reasonable
%    estimates of averaged posterior quantities. For example, do the
%    following to compute posterior inclusion probabilities (PIPs) averaged
%    over the hyperparameter settings:
%
%        w   = normalizelogweights(fit.logw);
%        PIP = fit.alpha * w(:);
%
%    And do the following to compute the posterior mean estimate of sa:
%
%        w       = normalizelogweights(fit.logw);
%        mean_sa = dot(fit.sa(:),w(:));
%
%    This is precisely how final posterior quantities are reported by
%    varbvsprint (type 'help varbvsprint' for more details). To account for
%    discrepancies between the prior on (sigma,sa,logodds) and the sampling
%    density used to draw candidate settings of the hyperparameters, adjust
%    the log-marginal probabilities weights by setting fit.logw = fit.logw +
%    logp/logq, where logp is the log-density of the prior distribution, and
%    logq is the log-density of the sampling distribution. (This is
%    importance sampling; see, for example, R. M. Neal, "Annealed importance
%    sampling", Statistics and Computing, 2001.)
%
% PRIOR ON PROPORTION OF VARIANCE EXPLAINED:
%    Specifying the prior variance of the regression coefficients (sa) can
%    be difficult, which is why we have included the option of fitting this
%    hyperparameter to the data (see input argument update_sa
%    above). However, in many settings, especially when a small number of
%    variables are included in the regression model, it is preferrable to
%    average over candidate settings of sa instead of fitting sa to the
%    data. To choose a set of candidate settings for sa, we have advocated
%    for setting sa indirectly through a prior estimate of the proportion of
%    variance in the outcome explained by the variables (abbreviated as
%    PVE), since it is often more natural to specify the PVE rather than the
%    prior variance (see references below). This is technically only
%    suitable or the linear regression model (family = 'gaussian'), but
%    could potentially be used for the linear regression model in an
%    approximate way.
%
%    For example, one could approximate a uniform prior on the PVE by
%    drawing the PVE uniformly between 0 and 1, additionally specifying
%    candidate settings for the prior log-odds, then computing the prior
%    variance (sa) as follows:
%
%      X  = X - repmat(mean(X),size(X,1),1);
%      sx = sum(var1(X));
%      sa = PVE./(1-PVE)./(sigmoid(log(10)*logodds)*sx)}
% 
%    Note that this calculation will yield sa = 0 when PVE = 0, and sa = Inf
%    when PVE = 1. The first line sets the mean of each column of X to 0,
%    which is needed to ensure that var1 correctly computes the sample
%    variance of each variable.
% 
%    Also, bear in mind that if there are additional covariates (Z) included
%    in the linear regression model that explain variance in Y, then it will
%    usually make more sense to first remove the linear effects of these
%    covariates before performing this calculation. The PVE would then
%    represent the prior proportion of variance in the residuals of Y that
%    are explained by the candidate variables. For an example of how to do
%    this, see varbvs.m, under "preprocessing steps". Alternatively, one
%    could include the matrix Z in the calculation above, taking care to
%    ensure that the covariates are included in the model with
%    probability 1.
%
% MEMORY REQUIREMENTS:
%    Finally, we point out that the optimization procedures were carefully
%    designed so that they can be applied to very large data sets; to date,
%    this code has been tested on data sets with >500,000 variables and
%    >10,000 samples. An important limiting factor is the ability to store
%    the data matrix X in memory. To reduce memory requirements, in the
%    MATLAB interface we require that X be single precision (type 'help
%    single'). Additionally, we mostly avoid generating intermediate
%    products that are of the same size as X. Only one such intermediate
%    product is generated when family = 'gaussian', and none for family =
%    'binomial'.
%
% LICENSE: GPL v3
%
% DATE: February 21, 2016
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
%    Y. Guan, M. Stephens (2011). Bayesian variable selection regression for
%    genome-wide association studies and other large-scale problems. Annals
%    of Applied Statistics 5, 1780–-1815.
% 
%    X. Zhou, P. Carbonetto, M. Stephens (2013). Polygenic modeling with
%    Bayesian sparse linear mixed models. PLoS Genetics 9, e1003264.
%
% SEE ALSO:
%    varbvsprint, varbvscoefcred
%
% EXAMPLES:
%    See demo_qtl.m and demo_cc.m for examples.
%
function fit = varbvs (X, Z, y, labels, family, options)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % (1) CHECK INPUTS
  % ----------------
  % Input X must be single precision, and cannot be sparse.
  if issparse(X)
    error('Input X cannot be sparse');
  end
  if ~isa(X,'single')
    X = single(X);
  end

  % If input Z is not empty, it must be double precision, and must have as
  % many rows as X.
  if ~isempty(Z)
    if size(Z,1) ~= n
      error('Inputs X and Z do not match.');
    end
    Z = double(full(Z));
  end

  % Add intercept.
  Z    = [ones(n,1) Z];
  ncov = size(Z,2) - 1;

  % Input y must be a double-precision column vector with n elements.
  y = double(y(:));
  if length(y) ~= n
    error('Inputs X and y do not match.');
  end

  % The labels must be a cell array with p elements, or an empty array.
  if nargin < 4
    labels = [];
  end
  if isempty(labels)
    labels = cellfun(@num2str,num2cell(1:p)','UniformOutput',false);
  else
    labels = labels(:);
    if (~iscell(labels) | length(labels) ~= p)
      error('labels must be a cell array with numel(labels) = size(X,2)');
    end
  end
  
  % By default, Y is a quantitative trait (normally distributed).
  if nargin < 5
    family = 'gaussian';
  end
  if isempty(family)
    family = 'gaussian';
  end
  if ~(strcmp(family,'gaussian') | strcmp(family,'binomial'))
    error('family must be gaussian or binomial');
  end
  
  % (2) PROCESS OPTIONS
  % -------------------
  % If the 'options' input argument is not specified, all the options are
  % set to the defaults.
  if nargin < 6
    options = [];
  end
  
  % OPTIONS.TOL
  % Set the convergence tolerance of the co-ordinate ascent updates.
  if isfield(options,'tol')
    tol = options.tol;
  else
    tol = 1e-4;
  end

  % OPTIONS.MAXITER
  % Set the maximum number of inner-loop iterations.
  if isfield(options,'maxiter')
    maxiter = options.maxiter;
  else
    maxiter = 1e4;
  end
  if ~isfinite(maxiter)
    error('options.maxiter must be a finite number');
  end

  % OPTIONS.VERBOSE
  % Determine whether to output progress to the console.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end

  % OPTIONS.SIGMA
  % Get candidate settings for the variance of the residual, if provided.
  % Note that this option is not valid for a binary trait.
  if isfield(options,'sigma')
    sigma        = double(options.sigma(:)');
    update_sigma = false;
    if strcmp(family,'binomial')
      error('options.sigma is not valid with family = binomial');
    end
  else
    sigma        = var(y);
    update_sigma = true;
  end

  % OPTIONS.SA
  % Get candidate settings for the prior variance of the coefficients, if
  % provided.
  if isfield(options,'sa')
    sa        = double(options.sa(:)');
    update_sa = false;
  else
    sa        = 1;
    update_sa = true;
  end

  % OPTIONS.LOGODDS
  % Get candidate settings for the prior log-odds of inclusion. This may
  % either be specified as a vector, in which case this is the prior applied
  % uniformly to all variables, or it is a p x ns matrix, where p is the
  % number of variables and ns is the number of candidate hyperparameter
  % settings, in which case the prior log-odds is specified separately for
  % each variable. A default setting is only available if the number of
  % other hyperparameter settings is 1, in which case we select 20 candidate
  % settings for the prior log-odds, evenly spaced between log10(1/p) and
  % -1.
  if isfield(options,'logodds')
    logodds = double(options.logodds);
  elseif isscalar(sigma) & isscalar(sa)
    logodds = linspace(-log10(p),-1,20);
  else
    error('options.logodds must be specified');
  end
  if size(logodds,1) == p
    prior_same = false;
  else
    prior_same = true;
    logodds    = logodds(:)';
  end

  % Here is where I ensure that the numbers of candidate hyperparameter
  % settings are consistent.
  ns = max([numel(sigma) numel(sa) size(logodds,2)]);
  if isscalar(sigma)
    sigma = repmat(sigma,1,ns);
  end
  if isscalar(sa)
    sa = repmat(sa,1,ns);
  end
  if size(logodds,2) == 1
    logodds = repmat(logodds,1,ns);
  end
  if numel(sigma) ~= ns | numel(sa) ~= ns | size(logodds,2) ~= ns
    error('options.sigma, options.sa and options.logodds are inconsistent');
  end

  % OPTIONS.UPDATE_SIGMA
  % Determine whether to update the residual variance parameter. Note
  % that this option is not valid for a binary trait.
  if isfield(options,'update_sigma')
    update_sigma = options.update_sigma;
    if strcmp(family,'binomial')
      error('options.update_sigma is not valid with family = binomial');
    end
  end

  % OPTIONS.UPDATE_SA
  % Determine whether to update the prior variance of the regression
  % coefficients.
  if isfield(options,'update_sa')
    update_sa = options.update_sa;
  end

  % OPTIONS.SA0
  % Get the scale parameter for the scaled inverse chi-square prior.
  if isfield(options,'sa0')
    sa0 = options.sa0;
  else
    sa0 = 0;
  end

  % OPTIONS.N0
  % Get the number of degrees of freedom for the scaled inverse chi-square
  % prior.
  if isfield(options,'n0')
    n0 = options.n0;
  else
    n0 = 0;
  end

  % OPTIONS.ALPHA
  % Set initial estimates of variational parameter alpha.
  initialize_params = true;
  if isfield(options,'alpha')
    alpha = double(options.alpha);
    initialize_params = false;  
    if size(alpha,1) ~= p
      error('options.alpha must have as many rows as X has columns');
    end
    if size(alpha,2) == 1
      alpha = repmat(alpha,1,ns);
    end
  else
    alpha = rand(p,ns);
    alpha = alpha ./ repmat(sum(alpha),p,1);
  end

  % OPTIONS.MU
  % Set initial estimates of variational parameter mu.
  if isfield(options,'mu')
    mu = double(options.mu);
    initialize_params = false;  
    if size(mu,1) ~= p
      error('options.mu must have as many rows as X has columns');
    end
    if size(mu,2) == 1
      mu = repmat(mu,1,ns);
    end
  else
    mu = randn(p,ns);
  end

  % OPTIONS.ETA
  % Set initial estimates of variational parameter eta. Note this is only
  % relevant for logistic regression.
  if isfield(options,'eta')
    eta               = double(options.eta);
    initialize_params = false;
    optimize_eta      = false;
    if ~strcmp(family,'binomial')
      error('options.eta is only valid for family = binomial');
    end
    if size(eta,1) ~= n
      error('options.eta must have as many rows as X');
    end
    if size(eta,2) == 1
      eta = repmat(eta,1,ns);
    end
  else
    eta          = ones(n,ns);
    optimize_eta = true;
  end

  % OPTIONS.OPTIMIZE_ETA
  % Determine whether to update the variational parameter eta. Note this
  % is only relevant for logistic regression.
  if isfield(options,'optimize_eta')
    optimize_eta = options.optimize_eta;
    if ~strcmp(family,'binomial')
      error('options.optimize_eta is only valid for family = binomial');
    end
  end
  
  % OPTIONS.INITIALIZE_PARAMS
  % Determine whether to find a good initialization for the variational
  % parameters.
  if isfield(options,'initialize_params')
    initialize_params = options.initialize_params;
    if initialize_params & ns == 1
      error(['options.initialize_params = true has no effect when ' ...
             'there is only one hyperparameter setting.']);
    end
  end

  % OPTIONS.NR
  % This is the number of samples to draw of the proportion of variance
  % in Y explained by the fitted model.
  if isfield(options,'nr')
    nr = options.nr;
  else
    nr = 100;
  end
  
  % TO DO: Allow specification of summary statistics ('Xb') from "fixed"
  % variational estimates for an external set of variables.
  clear options

  % (3) PREPROCESSING STEPS
  % -----------------------
  % Adjust the genotypes and phenotypes so that the linear effects of
  % the covariates are removed. This is equivalent to integrating out
  % the regression coefficients corresponding to the covariates with
  % respect to an improper, uniform prior; see Chipman, George and
  % McCulloch, "The Practical Implementation of Bayesian Model
  % Selection," 2001.
  if strcmp(family,'gaussian')

    % Here I compute two quantities that are used here to remove linear
    % effects of the covariates (Z) on X and y, and later on (in function 
    % "outerloop"), to efficiently compute estimates of the regression
    % coefficients for the covariates.
    SZy = (Z'*Z)\(Z'*y);
    SZX = double((Z'*Z)\(Z'*X));
    
    if ncov == 0
      X = X - repmat(mean(X),length(y),1);
      y = y - mean(y);
    else

      % This should give the same result as centering the columns of X and
      % subtracting the mean from y when we have only one covariate, the
      % intercept.
      y = y - Z*SZy;
      X = X - Z*SZX;
    end
  else
    SZy = [];
    SZX = [];
  end

  % Provide a brief summary of the analysis.
  if verbose
    fprintf('Welcome to           ');
    fprintf('--       *                              *               \n');
    fprintf('VARBVS version 2.0.3 ');
    fprintf('--       |              |               |               \n');
    fprintf('large-scale Bayesian ');
    fprintf('--       ||           | |    |          || |     |   |  \n');
    fprintf('variable selection   ');
    fprintf('-- |     || | |    |  | ||  ||        |||| ||    |   || \n');
    fprintf('*********************');
    fprintf('********************************************************\n');
    fprintf('Copyright (C) 2012-2016 Peter Carbonetto.\n')
    fprintf('See http://www.gnu.org/licenses/gpl.html for the full ');
    fprintf('license.\n');
    fprintf('Fitting variational approximation for Bayesian variable ');
    fprintf('selection model.\n');
    fprintf('family:     %-8s',family); 
    fprintf('   num. hyperparameter settings: %d\n',numel(sa));
    fprintf('samples:    %-6d',n); 
    fprintf('     convergence tolerance         %0.1e\n',tol);
    fprintf('variables:  %-6d',p); 
    fprintf('     iid variable selection prior: %s\n',tf2yn(prior_same));
    fprintf('covariates: %-6d',ncov);
    fprintf('     fit prior var. of coefs (sa): %s\n',tf2yn(update_sa));
    fprintf('intercept:  yes        ');
    if strcmp(family,'gaussian')
      fprintf('fit residual var. (sigma):    %s\n',tf2yn(update_sigma));
    elseif strcmp(family,'binomial')
      fprintf('fit approx. factors (eta):    %s\n',tf2yn(optimize_eta));
    end
  end
  
  % (4) INITIALIZE STORAGE FOR THE OUTPUTS
  % --------------------------------------
  % Initialize storage for the variational estimate of the marginal log-
  % likelihood for each hyperparameter setting (logw), the variances of the
  % regression coefficients (s), and the posterior mean estimates of the
  % coefficients for the covariates, including the intercept (muz).
  logw   = zeros(1,ns);
  s      = zeros(p,ns);
  mu_cov = zeros(ncov+1,ns);

  % (5) FIT BAYESIAN VARIABLE SELECTION MODEL TO DATA
  % -------------------------------------------------
  if ns == 1

    % Find a set of parameters that locally minimize the K-L
    % divergence between the approximating distribution and the exact
    % posterior.
    if verbose
      fprintf('        variational    max.   incl variance params\n');
      fprintf(' iter   lower bound  change   vars   sigma      sa\n');
    end
    [logw sigma sa alpha mu s eta mu_cov] = ...
        outerloop(X,Z,y,family,SZy,SZX,sigma,sa,logodds,alpha,mu,eta,tol,...
                  maxiter,verbose,[],update_sigma,update_sa,optimize_eta,...
                  n0,sa0);
    if verbose
      fprintf('\n');
    end
  else
      
    % If a good initialization isn't already provided, find a good
    % initialization for the variational parameters. Repeat for each
    % candidate setting of the hyperparameters.
    if initialize_params
      if verbose
        fprintf('Finding best initialization for %d combinations of ',ns);
        fprintf('hyperparameters.\n');
        fprintf('-iteration-   variational    max.   incl variance params\n');
        fprintf('outer inner   lower bound  change   vars   sigma      sa\n');
      end
      
      % Repeat for each setting of the hyperparameters.
      for i = 1:ns
        [logw(i) sigma(i) sa(i) alpha(:,i) mu(:,i) s(:,i) eta(:,i) ...
         mu_cov(:,i)] = ...
            outerloop(X,Z,y,family,SZy,SZX,sigma(i),sa(i),logodds(:,i),...
                      alpha(:,i),mu(:,i),eta(:,i),tol,maxiter,verbose,i,...
                      update_sigma,update_sa,optimize_eta,n0,sa0);
      end
      if verbose
        fprintf('\n');
      end

      % Choose an initialization common to all the runs of the coordinate
      % ascent algorithm. This is chosen from the hyperparameters with
      % the highest variational estimate of the marginal likelihood.
      [ans i] = max(logw);
      alpha   = repmat(alpha(:,i),1,ns);
      mu      = repmat(mu(:,i),1,ns);
      if optimize_eta
        eta = repmat(eta(:,i),1,ns);
      end
      if update_sigma
        sigma = repmat(sigma(i),1,ns);
      end
      if update_sa
        sa = repmat(sa(i),1,ns);
      end
    end
    
    % Compute a variational approximation to the posterior distribution
    % for each candidate setting of the hyperparameters.
    if verbose
      fprintf('Computing marginal likelihood for %d combinations of ',ns);
      fprintf('hyperparameters.\n');
      fprintf('-iteration-   variational    max.   incl variance params\n');
      fprintf('outer inner   lower bound  change   vars   sigma      sa\n');
    end

    % For each setting of the hyperparameters, find a set of parameters that
    % locally minimize the K-L divergence between the approximating
    % distribution and the exact posterior.
    for i = 1:ns
      [logw(i) sigma(i) sa(i) alpha(:,i) mu(:,i) s(:,i) eta(:,i) ...
       mu_cov(:,i)]=...
          outerloop(X,Z,y,family,SZy,SZX,sigma(i),sa(i),logodds(:,i),...
                    alpha(:,i),mu(:,i),eta(:,i),tol,maxiter,verbose,i,...
                    update_sigma,update_sa,optimize_eta,n0,sa0);
    end
    if verbose
      fprintf('\n');
    end
  end

  % (6) CREATE FINAL OUTPUT
  % -----------------------
  if strcmp(family,'gaussian')
    fit = struct('family',family,'n',n,'labels',{labels},'n0',n0,'sa0',sa0,...
                 'mu_cov',{mu_cov},'update_sigma',update_sigma,'update_sa',...
                 update_sa,'prior_same',prior_same,'logw',{logw},...
                 'sigma',{sigma},'sa',sa,'logodds',{logodds},'alpha',...
                 {alpha},'mu',{mu},'s',{s},'eta',[],'optimize_eta',false);

    % Compute the proportion of variance in Y, after removing linear
    % effects of covariates, explained by the regression model.
    if verbose
      fprintf('Estimating proportion of variance in Y explained by model.\n');
    end
    fit.model_pve = varbvspve(X,fit,nr);

    % Compute the proportion of variance in Y, after removing linear
    % effects of covariates, explained by each variable.
    fit.pve = zeros(p,ns);
    sx      = var1(X);
    for i = 1:ns
      fit.pve(:,i) = sx.*(mu(:,i).^2 + s(:,i))/var(y,1);
    end
  elseif strcmp(family,'binomial')
    fit = struct('family',family,'n',n,'labels',{labels},'n0',n0,'sa0',sa0,...
                 'mu_cov',{mu_cov},'update_sa',update_sa,'optimize_eta',...
                 optimize_eta,'prior_same',prior_same,'logw',{logw},...
                 'sa',sa,'logodds',{logodds},'alpha',{alpha},'mu',{mu},...
                 's',{s},'eta',{eta},'update_sigma',false,'pve',[]);
  end
  
% ------------------------------------------------------------------
% This function implements one iteration of the "outer loop".
function [logw, sigma, sa, alpha, mu, s, eta, mu_cov] = ...
        outerloop (X, Z, y, family, SZy, SZX, sigma, sa, logodds, alpha, ...
                   mu, eta, tol, maxiter, verbose, outer_iter, ...
                   update_sigma, update_sa, optimize_eta, n0, sa0)
  p = length(alpha);
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  end

  % Note that we need to multiply the prior log-odds by log(10), because
  % varbvsnorm, varbvsbin and varbvsbinz define the prior log-odds using
  % the natural logarithm (base e).
  if strcmp(family,'gaussian')

    % Optimize the variational lower bound for the Bayesian variable
    % selection model.
    [logw err sigma sa alpha mu s] = ...
        varbvsnorm(X,y,sigma,sa,log(10)*logodds,alpha,mu,tol,maxiter,...
                   verbose,outer_iter,update_sigma,update_sa,n0,sa0);

    % Adjust the variational lower bound to account for integral over the
    % regression coefficients corresponding to the covariates.
    logw = logw - logdet(Z'*Z)/2;

    % Compute the posterior mean estimate of the regression
    % coefficients for the covariates under the current variational 
    % approximation.
    mu_cov = SZy - SZX*(alpha.*mu);
  elseif strcmp(family,'binomial')

    % Optimize the variational lower bound for the Bayesian variable
    % selection model.
    if size(Z,2) == 1
      [logw err sa alpha mu s eta] = ...
          varbvsbin(X,y,sa,log(10)*logodds,alpha,mu,eta,tol,maxiter,verbose,...
                    outer_iter,update_sa,optimize_eta,n0,sa0);
    else
        [logw err sa alpha mu s eta] = ...
            varbvsbinz(X,Z,y,sa,log(10)*logodds,alpha,mu,eta,tol,maxiter,...
                       verbose,outer_iter,update_sa,optimize_eta,n0,sa0);
    end

    % Compute the posterior mean estimate of the regression
    % coefficients for the covariates under the current variational
    % approximation.
    Xr     = X*(alpha.*mu);
    d      = slope(eta);
    mu_cov = (Z'*diag(sparse(d))*Z)\(Z'*(y - 0.5 - d.*Xr));
  end
  logw = logw(end);

% ------------------------------------------------------------------
function y = tf2yn (x)
  if x
    y = 'yes';
  else
    y = 'no';
  end
