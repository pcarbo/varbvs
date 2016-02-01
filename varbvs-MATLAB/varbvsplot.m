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
%          options.groups
%          options.gap
%          options.vars
%          options.pip
%
% DETAILS:
%    <Details go here>
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
%    varbvs.
%

% TO DO: Add description for this function here (see varbvsprint.m).
%
% NOTE: Order of variables shown is based on order in which the groups
% are assigned to the variables.
%
function varbvsplot (fit, options)

  % Some plot settings.
  marker_color = [0.10 0.10 0.44];  % Midnight blue.
  label_color  = [1.00 0.00 1.00];  % Magenta.
  marker_size  = 6;
    
  % PROCESS OPTIONS
  % ---------------
  % If the 'options' input argument is not specified, all the options are set
  % to the defaults.
  if nargin < 2
    options = [];
  end

  % OPTIONS.PIP
  % Calculate the posterior inclusion probabilities (PIPs) if they aren't
  % provided as one of the inputs.
  if isfield(options,'pip')
    pip = options.pip;
  else
    w   = normalizelogweights(fit.logw);
    pip = fit.alpha * w(:);
  end
  p = numel(pip);
  
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
    groups       = options.groups;
    group_labels = unique(groups,'stable');
    group_labels = group_labels(:)';
  else
    groups = ones(p,1);
  end

  % OPTIONS.GAP
  % Determine how much whitespace to draw in between each group of
  % variables. 
  if isfield(options,'gap')
    gap = options.gap;
  else
    gap = 0;
  end
  clear options

  % (3) GENERATE GENOME-WIDE SCAN PLOT
  % ----------------------------------
  % Draw the posterior probabilities. Repeat for each of the groups.
  pos    = 0;
  xticks = [];
  for i = group_labels

    % Plot the variables assigned to the ith group.
    j = find(groups == i);
    m = length(j);
    plot(pos + (1:m),pip(j),'o','MarkerFaceColor',marker_color,...
         'MarkerEdgeColor','none','MarkerSize',marker_size);

    % Highlight the selected variables, and add labels to them.
    hold on
    [j k] = intersect(j,vars);
    plot(pos + k,pip(j),'o','MarkerFaceColor',label_color,...
         'MarkerEdgeColor','none','MarkerSize',marker_size);
    text(pos + k,pip(j),strcat(repmat({' '},length(j),1),fit.labels(j)),...
         'Color',label_color,'HorizontalAlignment','left',...
         'VerticalAlignment','bottom','FontSize',10);

    % Add the group label.
    xticks = [xticks pos+m/2];
    
    pos = pos + m + gap;
  end
  hold off
  a = gca;
  set(a,'XLim',[0 pos-gap+1],'XTick',xticks,'XTickLabel',group_labels);
  set(a,'TickLength',[0.005 0.005],'FontSize',12);
  set(a,'TickDir','out');
  drawnow;
  xruler = a.XRuler;
  yruler = a.YRuler;
  xruler.Axle.Visible = 'off';
  yruler.Axle.Visible = 'off';
