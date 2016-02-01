% TO DO: Add description for this function here (see varbvsprint.m).
function varbvsplot (fit, options)

  % Plot settings.
  marker_color = [0.10 0.10 0.44];  % Midnight blue.
  label_color  = [1.00 0.55 0.00];  % Dark orange.
  marker_size  = 6;
    
  % PROCESS OPTIONS
  % ---------------
  % If the 'options' input argument is not specified, all the options are set
  % to the defaults.
  if nargin < 2
    options = [];
  end

  % OPTIONS.PIP
  % -----------
  if isfield(options,'pip')
    pip = options.pip;
  else
    w   = normalizelogweights(fit.logw);
    pip = fit.alpha * w(:);
  end
  p = numel(pip);
  
  % OPTIONS.N
  % The number of variables to label in the plot.
  if isfield(options,'n')
    n = options.n;
  else
    n = 0;
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
  % Get the variables to be labeled in the plot.
  if n == 0
    vars = [];
  else
    [ans vars] = sort(-pip);
    vars       = vars(1:n);
    vars       = vars(:)';
  end
  
  % Draw the posterior probabilities. Repeat for each group.
  pos = 0;
  for i = group_labels

    % Plot the variables assigned to group i.
    j = find(groups == i);
    m = length(j);
    plot(pos + (1:m),pip(j),'o','MarkerFaceColor',marker_color,...
         'MarkerEdgeColor','none','MarkerSize',marker_size);

    % Add the group label.
    text(pos + m/2,-0.125,num2str(i),'Color','black','HorizontalAlignment',...
         'center','VerticalAlignment','bottom','FontSize',12);

    % Add labels to the variables.
    [j k] = intersect(j,vars);
    text(pos + k,pip(j),strcat(repmat({' '},length(j),1),fit.labels(j)),...
         'Color',label_color,'HorizontalAlignment','left',...
         'VerticalAlignment','bottom','FontSize',10);
    
    pos = pos + m + gap;
    hold on
  end
  hold off
  a      = gca;
  set(a,'TickDir','out','YLim',[-0.15 1.05],'YTick',0:0.25:1,'XTick',[]);
  set(a,'XLim',[0 pos-gap+1]);
  set(a,'FontSize',12,'FontName','Arial');
  ylabel('posterior probability');
  drawnow;
  xruler = a.XRuler;
  yruler = a.YRuler;
  xruler.Axle.Visible = 'off';
  yruler.Axle.Visible = 'off';
  
  