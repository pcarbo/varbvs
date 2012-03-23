% In this small example, we explore the posterior distribution of the
% coefficients in a linear regression, with spike and slab priors. In this
% idealized case, the variables are independent. We compare two approaches
% to computing the posterior probabilities: one using a fully-factorized
% variational approximation, and a Markov chain Monte Carlo (MCMC) method.
% This script effectively runs a single trial for the first experiment
% described in the Bayesian Analysis paper.
clear

% SCRIPT PARAMETERS.
p  = 1e3;  % The number of variables (SNPs).
n  = 500;  % The number of samples.
na = 20;   % Number of variables that effect the outcome ("causal" SNPs).
ns = 2e4;  % Length of Markov chain.
m0 = 100;  % Upper limit on number of selected variables in Markov chain.

% The level of residual noise is set so that the proportion of variance
% explained should be a little less than 30% on average (according to my
% rough calculations).
se = 9;

% These two parameters specify the Beta prior on the proportion of variables
% (SNPs) that are included in the linear model of Y.
a = 0.02;
b = 1;

% This parameter specifies the prior on the variance of the regression
% coefficients (sa). For more information on this prior, see VARSIMBVS
% and BVSMCMC.
c = 0.02;

% Candidate values of the variance of the residual (sigma), the prior
% variance of the regression coefficients (sa), and logarithm of the prior
% inclusion probability (log10q).
sigma  = (8:13)';
sa     = (0.025:0.025:0.4)';
log10q = (-2.5:0.25:-1)';

% Set the random number generator seed.
seed = 1;
rand('state',seed);
randn('state',seed);

% CREATE THE DATA.
% Note that X and y are centered.
fprintf('Creating data.\n');
[maf beta] = createsnps(p,na);
[X y]      = createdata(maf,beta,se,n);

% Calculate the proportion of variance explained. Here, SZ is the sample
% genetic variance.
sz = var(X*beta,1);
fprintf('Proportion of variance explained is %0.3f.\n',sz/(sz + se));

% COMPUTE VARIATIONAL ESTIMATES.
fprintf('Computing variational estimates.\n');
[sigma sa log10q] = ndgrid(sigma,sa,log10q);
[w alpha mu] = varsimbvs(X,y,sigma,sa,log10q,a,b,c);

% COMPUTE MCMC ESTIMATES.
fprintf('Computing MCMC estimates: ');
sx = sum(var1(X));
[ss sas qs PIP] = bvsmcmc(X,y,a,b,@(x) logpve(c*sx,x),m0,ns);
fprintf('\n');

% (4.) SHOW ERRORS IN MEAN ESTIMATES.
fprintf('                 Var    MCMC   diff.\n');

% Show posterior mean estimates of the residual variance (sigma).
x = dot(w(:),log10(sigma(:)));
y = mean(log10(ss));
fprintf('log10(sigma) %7.4f %7.4f %+0.4f\n',x,mean(y),x-y);

% Show posterior mean estimates of the prior variance (sa).
x = dot(w(:),log10(sa(:)));
y = mean(log10(sas));
fprintf('log10(sa)    %7.4f %7.4f %+0.4f\n',x,mean(y),x-y);

% Show posterior mean estimates of the prior inclusion probability (q).
x = dot(w(:),log10q(:));
y = mean(log10(qs));
fprintf('log10(q)     %7.4f %7.4f %+0.4f\n',x,mean(y),x-y);

% Show estimates of the correlation between the prior variance (sa) and
% the prior inclusion probability (q).
x = corrcoefw(log10(sa),log10q,w);
y = corrcoef(log10(sas),log10(qs));
y = y(2);
fprintf('                                               Var   MCMC\n');
fprintf('Correlation between log10(sa) and log10(q): %+0.3f %+0.3f\n',x,y);

% (5.) DISPLAY OTHER RESULTS.
clf
pos = get(gcf,'Position');
set(gcf,'Name','Simulation results','Color','white',...
	'NumberTitle','off','Position',[pos(1:2) 623 529]);

% Plot the variational estimate of the residual variance (sigma).
x = ndop(@mean,sigma,1);
N = ndsum(w,1);
subplot(3,3,1);
bar(x,N,1,'LineStyle','none','FaceColor',rgb('yellowgreen'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',[7 14],'XTick',8:2:12);
set(gca,'YLim',[0 0.8],'YTick',[0 0.4 0.8]);
set(gca,'TickDir','out','TickLength',[0.03 0.03]);
xlabel('\sigma^2');
ylabel('posterior');
title('Variational');

% Plot the MCMC estimate of the residual variance (sigma).
dx    = x(2) - x(1);
edges = [ x(1)-dx/2 ; x+dx/2];
N     = histc(ss,edges)/ns;
N     = N(1:end-1);
subplot(3,3,2);
bar(x,N,1,'LineStyle','none','FaceColor',rgb('dodgerblue'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',[7 14],'XTick',8:2:12);
set(gca,'YLim',[0 0.8],'YTick',[0 0.4 0.8]);
set(gca,'TickDir','out','TickLength',[0.03 0.03]);
xlabel('\sigma^2');
ylabel('posterior');
title('MCMC');

% Plot the variational estimate of the prior variance (sa).
x = ndop(@mean,sa,2);
N = ndsum(w,2);
x = (x(1:2:end) + x(2:2:end))/2;
N = N(1:2:end) + N(2:2:end);
subplot(3,3,4);
bar(x,N,1,'LineStyle','none','FaceColor',rgb('yellowgreen'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',[0 0.55]);
set(gca,'YLim',[0 0.8],'YTick',[0 0.4 0.8]);
set(gca,'TickDir','out','TickLength',[0.03 0.03]);
xlabel('\sigma_{\beta}^2');
ylabel('posterior');
title('Variational');

% Plot the MCMC estimate of the prior variance (sa).
dx    = x(2) - x(1);
edges = [ x(1)-dx/2; x+dx/2];
N     = histc(sas,edges)/ns;
N     = N(1:end-1);
subplot(3,3,5);
bar(x,N,1,'LineStyle','none','FaceColor',rgb('dodgerblue'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',[0 0.55]);
set(gca,'YLim',[0 0.8],'YTick',[0 0.4 0.8]);
set(gca,'TickDir','out','TickLength',[0.03 0.03]);
xlabel('\sigma_{\beta}^2');
ylabel('posterior');
title('MCMC');

% Plot the variational estimate of the prior inclusion probability (q).
x = ndop(@mean,log10q,3);
N = ndsum(w,3);
subplot(3,3,7);
bar(x,N,1,'LineStyle','none','FaceColor',rgb('yellowgreen'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',[-3 -0.5]);
set(gca,'YLim',[0 0.8],'YTick',[0 0.4 0.8]);
set(gca,'TickDir','out','TickLength',[0.03 0.03]);
xlabel('log_{10}\pi');
ylabel('posterior');
title('Variational');

% Plot the MCMC estimate of the prior inclusion probability (q).
dx    = x(2) - x(1);
edges = [ x(1)-dx/2 ; x+dx/2];
N     = histc(log10(qs),edges)/ns;
N     = N(1:end-1);
subplot(3,3,8);
bar(x,N,1,'LineStyle','none','FaceColor',rgb('dodgerblue'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',[-3 -0.5]);
set(gca,'YLim',[0 0.8],'YTick',[0 0.4 0.8]);
set(gca,'TickDir','out','TickLength',[0.03 0.03]);
xlabel('log_{10}\pi');
ylabel('posterior');
title('MCMC');

% Plot the ROC curve from the variational method.
subplot(3,3,6);
[fp tp] = createroc(beta ~= 0,alpha);
plot(fp,tp,'-','LineWidth',2,'Color',rgb('yellowgreen'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',na * [ -0.1 1.1 ]);
set(gca,'YLim',[ 0 max(tp)]);
xlabel('no. false +ves (\alpha)');
ylabel('no. true +ves (\beta)');
title('Variational');
box off

% Plot the ROC curve from the MCMC method.
subplot(3,3,9);
[fp tp] = createroc(beta ~= 0,PIP);
plot(fp,tp,'-','LineWidth',2,'Color',rgb('dodgerblue'));
set(gca,'FontSize',10,'FontName','Geneva');
set(gca,'XLim',na * [ -0.1 1.1 ]);
set(gca,'YLim',[ 0 max(tp)]);
xlabel('no. false +ves (\alpha)');
ylabel('no. true +ves (\beta)');
title('MCMC');
box off

