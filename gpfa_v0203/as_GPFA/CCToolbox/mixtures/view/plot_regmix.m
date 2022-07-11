function Options = plot_regmix(M,trajs,Options)
%Plot_Regmix  Plot regression mixture model
%
%   AxesHandles = Plot_Regmix(Model,[trajs],[Options])
%
%   If trajs is not empty, then the trajectory data will be plotted
%   using plot_seqdata(x, Y(:,d), ...) for each dimension d. The x
%   and Y come from [Y,x,Seq] = trajs2seq(trajs).
%
%   Options  (defaults in parentheses)
%   ----------------------
%   PlotMean - {0,1}, plot cluster means? (yes)
%   TwoD     - {0,1}, make a 2D plot as well as a 1D plot (yes)
%   delta    - delta value passed in to AXISCLOSE (1, pass in [] to bypass this)
%   PlotDimensions - {0,1}, plot separate dimensions? (yes)
%
% Comments:
%  - uses plot_regcoef; passes Options through.
%  - uses plot_seqdata; passes Options through.
%  - if not already set, plot_regcoef:len is set to 1+maxlen of trajs.


% Scott Gaffney   5 February 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'plot_regmix';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Begin Argument Processing
%%
trajs = cexist('trajs',[]);
Options = cexist('Options',[]);
M = SetFieldDef(M,'zero','nozero');
M.Options = SetFieldDef(M.Options,'MinLen',[]);
Options = SetFieldDef(Options,'PlotAlign',1);
Options = SetFieldDef(Options,'PlotMean',1);
Options = SetFieldDef(Options,'PlotSymbols',0);
Options = SetFieldDef(Options,'PlotDimensions',1);
Options = SetFieldDef(Options,'LineWidth',3);
Options = SetFieldDef(Options,'TwoD',0);
Options = SetFieldDef(Options,'delta',[5 10]);
Options = SetFieldDef(Options,'ScaleInSpace',0);
Options = SetFieldDef(Options,'AddPlot',0);
Options = SetFieldDef(Options,'SeparateFigs',0);
Options = SetFieldDef(Options,'PlotTrajsDayMarkers',0);
Options = SetFieldDef(Options,'PlotMeansDayMarkers',0);
Options = SetFieldDef(Options,'MeanCurveColor',[]);
Options = SetFieldDef(Options,'plot_regcoef',[]);

if (~isempty(trajs))
  [Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
  if (Options.PlotAlign)
    [x,Y] = AlignMat(M,x,Seq,Y);
    if (Options.ScaleInSpace)
      X = BasisMat(M,x);
    end
  end
  % set adaptive domain if domain is not prespecified
  if (~isfield(Options.plot_regcoef,'T'))
    Options.plot_regcoef.T = linspace(min(min(x)),max(max(x)),20);
  end
else               
  Y=[];  % tag for not plotting data
end
%%
%%% End Argument Processing



[P,K,D] = size(M.Mu);
if (~isempty(Y)), xd = size(x,2); end
child = sort(get(get(0,'CurrentFigure'),'children'));
lenc = length(child);

%%%%%
%%% Plot 2-D graph
if (Options.TwoD & D==2)
  if (Options.SeparateFigs)
    ah_2d = figure;  ah_2d = gca;
  else
    if (Options.PlotDimensions),  SubPlotDim = D+1;  else SubPlotDim = 1; end
    if (lenc~=SubPlotDim)
      ah_2d = subplot(1,SubPlotDim,SubPlotDim);
    else
      ah_2d = child(1);
      axes(ah_2d);
    end
    if (Options.AddPlot==0), cla; end
  end
  
  % plot the data
  if (~isempty(Y))
    if (Options.PlotTrajsDayMarkers), Options.PlotDayMarkers=1; end
    if (Options.ScaleInSpace & Options.PlotAlign)
      if (isfield(M,'Ee')), e=M.Ee;  else e=permute(M.e,[1 3 2]); end
      p = plot_seqdata_aa([],X,Y(:,1:2),Seq,M.Mu(:,:,1:2),M.Ee,M.C,Options);
    else
      p = plot_seqdata(Y(:,1),Y(:,2),Seq,M.C,Options);
    end
    Options.PlotDayMarkers = 0;
  end
  % plot mean coefficients
  if (Options.PlotMean)
    if (Options.PlotMeansDayMarkers), Options.PlotDayMarkers=1; end
    syms=Options.PlotSymbols;  Options.PlotSymbols=0;
    p = cell2vector(plot_regcoef(M.Mu,[],[],Options));
    set(p,'LineWidth',Options.LineWidth,'Tag','Means');
    if (~isempty(Options.MeanCurveColor))
      set(p,'color',Options.MeanCurveColor);
    end
    Options.PlotSymbols=syms;
    Options.PlotDayMarkers = 0;
  end
  SetAxesPos(ah_2d, [2 4], [.17 .75]);
  axisclose(ah_2d,Options.delta);
%   axisequal(ah_2d);

else  % otherwise set number of plots to the dimension Mu
  SubPlotDim = D;
  ah_2d = [];  % this is used below
end


%%%%%
%%% Plot the individual dimensions
ah=[];
if (Options.PlotDimensions)
  for d=1:D
    if (Options.SeparateFigs)
      fig = figure;   ah(end+1) = gca;
    else
      if ((Options.TwoD==0 & lenc~=D) | (Options.TwoD & lenc~=D+1))
        ah(end+1) = subplot(1,SubPlotDim,d);
      else
        if (Options.TwoD) 
          ah(end+1) = child(d+1);
        else
          ah(end+1) = child(d);
        end
        axes(ah(end));
      end
      if (Options.AddPlot==0), cla; end
    end
    
    % plot the data
    if (~isempty(Y))
      if (Options.PlotTrajsDayMarkers), Options.PlotDayMarkers=1; end
      shd = min(xd,d);
      if (Options.ScaleInSpace & Options.PlotAlign)
        if (isfield(M,'Ee')), e=M.Ee;  else e=permute(M.e,[1 3 2]); end
        shareE = size(e,3);
        she = min(shareE,d);
        p = plot_seqdata_aa(x(:,shd),X(:,:,shd),Y(:,d),Seq, ...
          M.Mu(:,:,d),e(:,:,she),M.C,Options);
      else
        p = plot_seqdata(x(:,shd),Y(:,d),Seq,M.C,Options);
      end
      Options.PlotDayMarkers = 0;
    end
    % plot the coefficients
    if (Options.PlotMean)
      if (Options.PlotMeansDayMarkers), Options.PlotDayMarkers=1; end
      syms=Options.PlotSymbols;  Options.PlotSymbols=0;
      p = cell2vector(plot_regcoef(M.Mu(:,:,d),[],[],Options));
      set(p,'LineWidth',Options.LineWidth,'Tag','Means');
      if (~isempty(Options.MeanCurveColor))
        set(p,'color',Options.MeanCurveColor);
      end
      Options.PlotSymbols=syms;
      Options.PlotDayMarkers = 0;
    end
  end
end

SetAxesPos(ah, [2 4], [.17 .75]);
Options.hnds = [ah_2d fliplr(ah)];

%SpreadAxes([ah ah_2d]);
%EqualLimits(ah,Options.delta);  % this puts the axes on the same scale
%axis(ah,'square');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BasisMat
%
function X = BasisMat(M,x)
[N,D] = size(x);
P = size(M.Mu,1);
X = zeros(N,P,D);

if (isfield(M,'knots'))
  for d=1:D
    X(:,:,d) = bsplinebasis(M.knots,M.order,x(:,d));
  end
else
  for d=1:D
    X(:,:,d) = regmat(x(:,d),P-1);
  end
end

