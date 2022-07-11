function Options = plot_gmix(M,trajs,Options)
%PLOT_GMIX  Plot a regression mixture model.
%
%   PLOT_GMIX(Model,[trajs],[Options])
%
%   If trajs is not empty, then the trajectory data will be plotted
%   using plot_seqdata(x, Y(:,d), ...) for each dimension d. The x and
%   Y come from [Y,x,Seq] = trajs2seq(trajs). If you do not want to plot
%   the data (i.e., trajs is empty), then the mean curves are plotted
%   at T = (0:m-1) where m = size(M.Mu,1). If you do not want this
%   behavior, then you can pass in a new value of T as Options.plot_gmix.T.
%
%   Options  (defaults in parentheses)
%   ----------------------
%   PlotMean - {0,1}, plot cluster means? (yes)
%   TwoD     - {0,1}, make a 2D plot as well as a 1D plot (yes)
%   plot_gmix.T - alternate vector containing plotting points for mean curves
%   delta    - delta value passed in to AXISCLOSE (1, pass in [] to bypass this)
%
%
% Comments:
%  - uses plot_regcoef; passes Options through.
%  - uses plot_seqdata; passes Options through.
%  - if not already set, plot_regcoef:len is set to maxlen of trajs.


% Scott Gaffney   5 February 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'plot_gmix';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Begin Argument Processing
%%
trajs = cexist('trajs',[]);
Options = cexist('Options',[]);
Options = SetFieldDef(Options,'PlotMean',1);
Options = SetFieldDef(Options,'TwoD',0);
Options = SetFieldDef(Options,'T',[]);
Options = SetFieldDef(Options,'delta',1);
Options = SetFieldDef(Options,'AddPlot',0);
Options = SetFieldDef(Options,'PlotTrajsDayMarkers',0);
Options = SetFieldDef(Options,'PlotMeansDayMarkers',0);
Options = SetFieldDef(Options,'MeanCurveColor',[]);
Options = SetFieldDef(Options,'PlotSymbols',0);
Options = SetFieldDef(Options,'LineWidth',3);
Options = SetFieldDef(Options,'SelectClass',[]);
Options = SetFieldDef(Options,'plot_gmix',[]);
Options.plot_gmix = SetFieldDef(Options.plot_gmix,'T',[]);
M = SetFieldDef(M,'Options',[]);
M.Options = SetFieldDef(M.Options,'MinLen',[]);
M = SetFieldDef(M,'zero','nozero');

ni = size(M.Mu,1);
if (~isempty(trajs))
  [Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
else
  Y=[];  % tag for not plotting data
  x = (0:ni-1)';
end
if (isfield(Options.plot_gmix,'T') & ~isempty(Options.plot_gmix.T))
  T = Options.plot_gmix.T(1:ni);
else
  T = x(1:ni);
end

%%
%%% End Argument Processing



[m,K,D] = size(M.Mu);
child = get(gcf,'children');
lenc = length(child);

%%%%%
%%% Plot 2-D graph
if (Options.TwoD & D==2)
  SubPlotDim = 3;
  if (lenc~=D+1)
    ah_2d = subplot(1,SubPlotDim,SubPlotDim);
  else
    ah_2d = child(end);  % the first subplot created, should be last
    axes(ah_2d);
  end
  if (Options.AddPlot==0), cla; end
  xlabel('y_1'); ylabel('y_2');
  
  % plot the data
  if (~isempty(Y))
    p = plot_seqdata(Y(:,1),Y(:,2),Seq,M.C,Options);
  end
  % plot mean vectors
  if (Options.PlotMean)
    if (Options.PlotMeansDayMarkers), Options.PlotDayMarkers=1; end
    syms=Options.PlotSymbols;  Options.PlotSymbols=0;
    for k=1:K
      if (isempty(Options.SelectClass) | any(Options.SelectClass==k))
        p = colorplot(M.Mu(:,k,1),M.Mu(:,k,2),k,Options);
        set(p,'LineWidth',Options.LineWidth,'Tag','Means');
        if (~isempty(Options.MeanCurveColor))
          set(p,'color',Options.MeanCurveColor);
        end
      end  
    end
    Options.PlotSymbols=syms;
    Options.PlotDayMarkers = 0;
  end
  SetAxesPos(ah_2d, [2 4], [.25 .65]);
  axisclose(ah_2d,Options.delta);
  axisequal(ah_2d);

else  % otherwise set number of plots to the dimension Mu
  SubPlotDim = D;
  ah_2d = [];  % this is used below
end


%%%%%
%%% Plot the individual dimensions
ah=[];
for d=1:D
  if ((Options.TwoD==0 & lenc~=D) | (Options.TwoD & lenc~=D+1))
    ah(end+1) = subplot(1,SubPlotDim,d);
  else
    ah(end+1) = child(d);
    axes(ah(end));
  end
  if (Options.AddPlot==0), cla; end
  if (d==1)  ylabel('y_1');  elseif (d==2)  ylabel('y_2'); end
  xlabel('Time');
  
  % plot the data
  if (~isempty(Y))
    p = plot_seqdata(x,Y(:,d),Seq,M.C,Options);
  end
  % plot the coefficients
  if (Options.PlotMean)
    if (Options.PlotMeansDayMarkers), Options.PlotDayMarkers=1; end
    syms=Options.PlotSymbols;  Options.PlotSymbols=0;
    for k=1:K
      if (isempty(Options.SelectClass) | any(Options.SelectClass==k))
        p = colorplot(T,M.Mu(:,k,d),k,Options);
        set(p,'LineWidth',Options.LineWidth,'Tag','Means');
        if (~isempty(Options.MeanCurveColor))
          set(p,'color',Options.MeanCurveColor);
        end
      end  
    end
    Options.PlotSymbols=syms;
    Options.PlotDayMarkers = 0;
  end
end

SetAxesPos(ah, [2 4], [.25 .65]);
Options.hnds = [ah_2d fliplr(ah)];
% SpreadAxes([ah ah_2d]);
% EqualLimits(ah,Options.delta);
%axis(ah,'square');
