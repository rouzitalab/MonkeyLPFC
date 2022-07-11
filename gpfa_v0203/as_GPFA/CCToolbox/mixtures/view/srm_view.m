function Options = srm_view(M,trajs,Options)
%SRM_View  View SRM_* model
%
%   AxesHandles = SRM_View(Model,[trajs],[Options])
%

% Scott Gaffney   5 February 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'srm_view';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

%%% Begin Argument Processing
%%
trajs = cexist('trajs',[]);
Options = cexist('Options',[]);
%%
%%% End Argument Processing

Options = DefaultOptions(Options);
bottompos = [0.1200    0.5676    0.3329    0.3267];
toppos    = [0.3900    0.5486    0.3743    0.3686];

% Activate the figure
if (~Options.SeparateFigs)
  if (isempty(Options.bfig))
    Options.bfig = figure('Units','normalized','Position',bottompos);
  end
  figure(Options.bfig);
end

Options.plot_regcoef.T = (M.knots(1):.1:M.knots(end))';
Options.plot_regcoef.X = bsplinebasis(M.knots,M.order,Options.plot_regcoef.T);
Options = plot_regmix(M,trajs,Options);

info = listmodels(M.method);
set(gcf,'Name',info.description);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%
function Options = DefaultOptions(Options);
Options = SetFieldDef(Options,'bfig',[]);
Options = SetFieldDef(Options,'SeparateFigs',0);
