%--------------------------------------------------------------------------
% varbvsplot.m: Summarize variable selection results in a single plot.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Generate a single plot that summarizes the results of fitting the
%    Bayesian variable selection model to the data. When the variables are
%    genetic markers and the groups are chromosomes, the figure resembles a
%    "Manhattan plot" typically used to summarize the results of a
%    genome-wide association study or quantitative trait locus (QTL)
%    mapping study.
%
% USAGE:
%    varbvsplot(fit)
%    varbvsplot(fit, options)
%
% INPUT ARGUMENTS:
% fit      Output of function varbvs.
% options  A structure (type 'help struct') specifying some plot
%          settings. More details about these options are given
%          below. Fields, with their default settings given, include:
%
%          options.groups (grouping of variables)
%          options.gap (size of gap between each group in plot)
%          options.vars (variables to highlight and label)
%          options.score (precomputed posterior probabilities, or "score")
%
% DETAILS:
%    Optional input options.groups specifies the grouping of the
%    variables. This must be an array with as many entry as variables, in
%    which each entry is a unique number specifying the group assignment. In
%    the plot, the groups are shown in the same order that they appear in
%    this array. By default, all variables are assigned to a single group.
%    options.gap specifies how much space to leave in between each group
%    of variables in the plot.
%
%    The variables are draw along the horizontal axis, grouped according to
%    options.groups. By default, the vertical axis shows the posterior
%    inclusion probability (PIP), averaged over the hyperparameter settings
%    (see 'help varbvsprint' for more details about how the averaged
%    posterior probabilities are calculated). An alternative statistic for
%    each variable may be provided for the vertical axis in options.score.
%    For example, to more closely reproduce a "Manhattan plot" for a
%    genome-wide association study, compute the PIPs that ignore
%    correlations between variables, then set options.score to these PIPs
%    on the log-scale:
%
%      options.score = log10(varbvsindep(fit,X,Z,y) * fit.w(:) + 1e-5);
%
%    Finally, options.vars can be used to highlight and label variables in
%    the plot. This must be an array containing the indices of the variables
%    to be highlighted. The labels used are those provided by fit.labels.
%    By default, no variables are highlighted.
%
% LICENSE: GPL v3
%
% DATE: February 1, 2016
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
% SEE ALSO:
%    varbvs, varbvsprint.
%
function varbvsplot (fit, options)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %

  % Some plot settings.
  marker_color = [0.10 0.10 0.44];  % Midnight blue.
  label_color  = [1.00 0.00 1.00];  % Magenta.
  marker_size  = 6;
    
  % PROCESS OPTIONS
  % ---------------
  % If the 'options' input argument is not specified, all the options are
  % set to the defaults.
  if nargin < 2
    options = [];
  end

  % OPTIONS.SCORE
  % Calculate the posterior inclusion probabilities (PIPs) if a "score"
  % isn't provided as one of the inputs. 
  if isfield(options,'score')
    y = options.score;
  else
    y = fit.pip;
  end
  p = numel(y);
  
  % OPTIONS.VARS
  % Get the variables to be highlighted and labeled in the plot.
  if isfield(options,'vars')
    vars = options.vars;
  else
    vars = [];
  end
  
  % OPTIONS.GROUP
  % Determine the grouping of the variables. By default, all the
  % variables are assigned to a single group.
  if isfield(options,'groups')
    groups = options.groups;
  else
    groups = ones(p,1);
  end
  group_labels = unique(groups,'stable');
  group_labels = group_labels(:)';

  % OPTIONS.GAP
  % Determine how much whitespace to draw in between each group of
  % variables. 
  if isfield(options,'gap')
    gap = options.gap;
  else
    gap = 0;
  end
  clear options

  % GENERATE GENOME-WIDE SCAN PLOT
  % ------------------------------
  % Determine the positions of the variables along the horizontal axis.
  x      = zeros(p,1);
  pos    = 0;
  xticks = [];
  for i = group_labels
    j      = find(groups == i);
    m      = length(j);
    x(j)   = pos + (1:m);
    xticks = [xticks pos+m/2];
    pos    = pos + m + gap;
  end

  % Plot the posterior probabilities.
  plot(x,y,'o','MarkerFaceColor',marker_color,'MarkerEdgeColor','none',...
       'MarkerSize',marker_size);

  % Highlight the selected variables, and add labels to them.
  if length(vars) > 0
    hold on
    plot(x(vars),y(vars),'o','MarkerFaceColor',label_color,...
         'MarkerEdgeColor','none','MarkerSize',marker_size);
    text(x(vars),y(vars),...
         strcat(repmat({' '},length(vars),1),fit.labels(vars)),...
         'Color',label_color,'HorizontalAlignment','left',...
         'VerticalAlignment','bottom','FontSize',10);
    hold off
  end
  set(gca,'XLim',[0 pos-gap+1],'XTick',xticks,'XTickLabel',group_labels);
  set(gca,'TickLength',[0.005 0.005],'FontSize',12);
  set(gca,'TickDir','out');

