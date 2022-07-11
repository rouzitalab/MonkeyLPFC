function [trajs,M] = simulate_trajs(n,method,M,len,Options)
%SIMULATE_TRAJS  Simulate from specified cluster model.
%   [TRAJS,M] = SIMULATE_TRAJS(N,M,METHOD,[LEN],[OPTIONS]) returns 
%   a cell array of length N of sampled trajectories of varying lengths (or
%   exactly of length LEN, if specified).
%
%   M is either a filename of the requested model to simulate 
%   from or the model itself. If M is a filename, it must be in the form
%   of one of the files in simulation/scripts.
%
%   METHOD must be one of the cluster models recognized by ListModels().
%
%   If Options.plotme==1, then the simulated model will be plotted. All
%   other options will be passed through.

% Scott J Gaffney   5 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'simulate_trajs';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

M = cexist('M',[]);
method = cexist('method',[]);
if (isempty(M) | isempty(method))
  error([PROGNAME,': you must supply a model and a method']);
end

trajs=[];
len = cexist('len',[]);
Options = cexist('Options',[]);
Options = SetFieldDef(Options,'plotme',0);
Options = SetFieldDef(Options,'bfig',[]);
Options = SetFieldDef(Options,'tfig',[]);
Options = SetFieldDef(Options,'RandGeneratorState',rand('state'));
Options = SetFieldDef(Options,'RandGeneratorNState',randn('state'));

% Try to resolve method
if (isstr(M))
  M = getmodel(M,PROGNAME);
  if (isempty(M)), return; end
end
M.method = method;

% Perform the sampling
save_state = rand('state');  save_nstate = randn('state');
rand('state',Options.RandGeneratorState);
randn('state',Options.RandGeneratorNState);
simfunc = listmodels(M.method,'simfun');
if (isempty(simfunc)), return; end
[trajs,M] = feval(simfunc,M,n,len,Options);

% Final stuff
M.zero = 'nozero';
M.C = M.TrueC;

% Plotting stuff
if (Options.plotme)
  Options = SetFieldDef(Options,'PlotData',1);
  Options = SetFieldDef(Options,'PlotSymbols',1);
  viewmodel(M,trajs,Options);
end

M.Options = Options;
rand('state',save_state');  randn('state',save_nstate);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function M = getmodel(fname,PROGNAME)
try M = eval(fname);
catch
  errorbox(lasterr,PROGNAME);
  M=[];
  return;
end  
