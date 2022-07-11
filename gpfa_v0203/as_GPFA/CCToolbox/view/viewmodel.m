function Options = viewmodel(model,trajs,Options)
%VIEWMODEL  Plot results of regression mixtures (lreg/hreg/shreg/etc.).
%
%   Options = VIEWMODEL(Model,[trajs],[Options])
%
%   Options
%   ----------------------------------
%   PlotCoef    - {0,1}, plot indivdiual regression coefficients
%   bfig        - handle to bottom-level figure
%
%
% Comments:
%  - sets plot_regmix:PlotMean to on
%  - uses plot_regmix; passes Options through.


% Scott Gaffney   05 April 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'viewmodel';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

%%% Begin Argument Processing
%%%
trajs = cexist('trajs',[]);
Options = cexist('Options',[]);
%%%
%%% End Argument Processing


% preprocessing
Options = DefaultOptions(Options);
model = SetFieldDef(model,'Options',[]);

% Handle random trajectory plotting
Options = CalcRandTraj(trajs,model.C,Options);

% Get specific 'view' method and call it
view_func = listmodels(model.method,'view');
Options = feval(view_func,model,trajs,Options);
drawnow;  % flush the event queue

% set up the 'handles' field
if (isempty(Options.hnds)), Options.hnds = gca; end
% set(Options.hnds,'box','on','FontSize',Options.FontSize);
% xhnds = get(Options.hnds,'Xlabel');
% yhnds = get(Options.hnds,'Ylabel');
% if (iscell(xhnds)), xhnds = cell2mat(xhnds); end
% if (iscell(yhnds)), yhnds = cell2mat(yhnds); end
% set(xhnds,'FontSize',Options.FontSize);
% set(yhnds,'FontSize',Options.FontSize);
%axis tight;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%
function Options = DefaultOptions(Options);
Options = SetFieldDef(Options,'PlotCoef',0);
Options = SetFieldDef(Options,'bfig',[]);
Options = SetFieldDef(Options,'hnds',[]);
% Options = SetFieldDef(Options,'tfig',[]);
Options = SetFieldDef(Options,'PlotMean',1);
Options = SetFieldDef(Options,'PlotNumRand',[]);
Options = SetFieldDef(Options,'SelectTrajs',[]);
Options = SetFieldDef(Options,'SeparateFigs',0);
Options = SetFieldDef(Options,'FontSize',14);
