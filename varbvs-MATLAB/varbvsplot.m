% TO DO: Add description for this function here (see varbvsprint.m).
%
% NOTE: Order of variables shown is based on order in which the groups
% are assigned to the variables.
%
function varbvsplot (fit, options)

  % Plot settings.
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
  % provided as input.
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
  % Determine the grouping of the variables.
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
  % Draw the posterior probabilities. Repeat for each group.
  pos    = 0;
  xticks = [];
  for i = group_labels

    % Plot the variables assigned to group i.
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
  set(a,'TickLength',[0.005 0.005],'FontSize',12,'FontName','Arial');
  set(a,'TickDir','out');
  drawnow;
  xruler = a.XRuler;
  yruler = a.YRuler;
  xruler.Axle.Visible = 'off';
  yruler.Axle.Visible = 'off';

  