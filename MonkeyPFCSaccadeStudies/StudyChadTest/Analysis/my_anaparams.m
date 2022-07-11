function analysisParams = my_anaparams(varargin)
%analysisParams = my_anaparams(analysis_name)
%
% analysis_name is an optional string argument to select a specific set of
% analysisParams, otherwise the default parameters are returned.
%
% analysisParams is a struct containing all the parameters necessary for a
% particular analysis, with the following fields:
% name              - Name of the set of analysis params. Searched when
%                     analysis_name optional argument is passed in.
% outcomeCodes      - Only trials with newOutcomeCodes that match one of
%                     the entries in this array will be included in the analysis.
% critFlip          - The name of the event after which to search for a saccade.
% binWidth          - Width (ms) of the time window used for binning when
%                     creating neural trajectories
% kernSD            - Width (ms) of the smoothing kernel used for neural trajectories
% classifyTarg      - Which trial outcome would you like to classify?
%                     Possible values: 'class' 'cueColour' 'sacPol' 'sacClass' 'targPol' 'targClass'
% minTrialsPerGroup - If one of your groups (e.g. all trials with sacClass == 1)
%                     has fewer than this many trials then exclude it.
% useUnits          - 'all' uses sorted and unsorted units. 'sorted' uses
%                     only unsorted units. 'merged' is akin to threshold
%                     crossing events.
% ptest             - 
% statPerm          - How many permutations to run when making a null distribution for stats.
% kcross            - number of cross validations
% anaWins           - an array of structs describing time windows in the
%                     trial to analyze
%
% Each anaWins entry has the following fields:
% name          - Name of this window. e.g. 'baseline' or 'delay'
% timeLockEvent - The trial event for t = 0. Possible values:
%                 fixationOnset targetOnset cueOnset cueOffset fixationOffset saccadeOnset
% winEdges      - The time window (in ms relative to timeLockEvent) to analyze.
% avoidEvent    - A trial event (with avoidWin) that must be excluded from this window.
% avoidWin      - The window (in ms) around avoidEvent that must be
%                 excluded from this analysis window.
% plotX         - When including this window in a plot, what should be its x-range?
% vBar          - When plotting, where to put vertical bars, if anywhere.
% colour        - When plotting, what colour should be ascribed to this window?

%% Specify the parameters of the different analysis windows.
% Baseline window.
ana_win = struct(...
    'name', 'baseline',...
    'timeLockEvent', 'targetOnset',...  % fixationOnset targetOnset cueOnset cueOffset fixationOffset saccadeOnset
    'winEdges', [-249 0],...
    'plotX', [-249 0],...
    'avoidEvent', 'cueOffset',...
    'avoidWin', [-100 Inf],...
    'vBar', 0,...
    'colour', 'k');

%Target response.
ana_win(end+1).name = 'visual';
ana_win(end).timeLockEvent = 'targetOnset';
ana_win(end).winEdges = [1 300];
ana_win(end).plotX = [1 300];
ana_win(end).avoidEvent = 'fixationOffset';
ana_win(end).avoidWin = [-200 Inf];
ana_win(end).vBar = 300;
ana_win(end).colour = 'b';

% Delay window.
ana_win(end+1).name = 'cue';
ana_win(end).timeLockEvent = 'targetOnset';
ana_win(end).winEdges = [301 1300];
ana_win(end).plotX = [301 1300];
ana_win(end).avoidEvent = 'fixationOffset';
ana_win(end).avoidWin = [-200 Inf];
ana_win(end).vBar = 1300;
ana_win(end).colour = 'm';

ana_win(end+1).name = 'delay';
ana_win(end).timeLockEvent = 'targetOnset';
ana_win(end).winEdges = [1301 1600];
ana_win(end).plotX = [1301 1600];
ana_win(end).avoidEvent = 'fixationOffset';
ana_win(end).avoidWin = [-50 Inf];
ana_win(end).vBar = nan;
ana_win(end).colour = 'g';

% Motor window.
ana_win(end+1).name = 'motor';
ana_win(end).timeLockEvent = 'saccadeOnset';
ana_win(end).winEdges = [-49 200];  % [-69 50] might be better for dPCA
ana_win(end).plotX = [1651 1900];
ana_win(end).avoidEvent = 'fixationOffset';
ana_win(end).avoidWin = [-Inf 0];
ana_win(end).vBar = 1700;
ana_win(end).colour = 'r';

% For MI time series
ana_win(end+1).name = 'timeSeriesTarget';
ana_win(end).timeLockEvent = 'targetOnset';
ana_win(end).winEdges = [-249 2400];
ana_win(end).avoidEvent = 'saccadeOnset';
ana_win(end).avoidWin = [Inf -Inf];
ana_win(end+1).name = 'timeSeriesSaccade';
ana_win(end).timeLockEvent = 'saccadeOnset';
ana_win(end).winEdges = [-2449 200];
ana_win(end).avoidEvent = 'targetOnset';
ana_win(end).avoidWin = [Inf -Inf];

% For classification
ana_win(end+1).name = 'trajectory';
ana_win(end).timeLockEvent = 'targetOnset';
ana_win(end).winEdges = [-249 1250];  % In this dataset, this is known to be present for all good trials.
ana_win(end).avoidEvent = 'saccadeOnset';
ana_win(end).avoidWin = [Inf -Inf];

%% Analysis Parameters
%Default
analysisParams = struct(...
    'name', 'default',...
    'outcomeCodes', 0,...
    'critFlip', 'fixationOffset',...
    'binWidth', 50,...
    'kernSD', 50,...
    'classifyTarg', 'sacPol',...  % See getNewClass.m; 'class' 'cueColour' 'sacPol' 'sacClass' 'targPol' 'targClass'
    'minTrialsPerGroup', 15,...
    'useUnits', 'sorted',...  % 'all' 'sorted' 'merged'
    'ptest', 0.1,...
    'statPerm', 1000,...  % Number of permutations for statistics
    'kcross', 10,...  % number of cross validations.
    'anaWins', ana_win(strcmpi({ana_win.name}, 'delay')));

%behaviour
analysisParams(end+1) = analysisParams(1);
analysisParams(end).name = 'behaviour';
analysisParams(end).outcomeCodes = [0 9];
analysisParams(end).critFlip = 'cueOffset';  % I'm OK with saccades that come after cueOffset but before fixationOffset
analysisParams(end).anaWins = [];

%tuning
analysisParams(end+1) = analysisParams(1);
analysisParams(end).name = 'tuning';
analysisParams(end).outcomeCodes = [0 9];
analysisParams(end).critFlip = 'cueOffset';
analysisParams(end).minTrialsPerGroup = 8;
[~, lib] = ismember({'baseline' 'visual' 'cue' 'delay' 'motor'}, {ana_win.name});
analysisParams(end).anaWins = ana_win(lib);

%timeseries
analysisParams(end+1) = analysisParams(strcmpi({analysisParams.name}, 'tuning'));
analysisParams(end).name = 'timeseries';
[~, lib] = ismember({'timeSeriesTarget' 'timeSeriesSaccade'}, {ana_win.name});
analysisParams(end).anaWins = ana_win(lib);

%dPCA
analysisParams(end+1) = analysisParams(strcmpi({analysisParams.name}, 'tuning'));
analysisParams(end).name = 'dPCA';
analysisParams(end).classifyTarg = 'combinedRegionRule';
analysisParams(end).outcomeCodes = 0;  % Only correct
% analysisParams(end).useUnits = 'all';
% analysisParams(end).binWidth = 20;
% analysisParams(end).kernSD = 30;
[~, lib] = ismember({'trajectory'}, {ana_win.name});
analysisParams(end).anaWins = ana_win(lib);

%For classification - uses almost all windows.
analysisParams(end+1) = analysisParams(strcmpi({analysisParams.name}, 'tuning'));
analysisParams(end).name = 'classification';
% analysisParams(end).useUnits = 'all';
[~, lib] = ismember({'baseline' 'visual' 'cue' 'delay' 'motor' 'trajectory'}, {ana_win.name});
analysisParams(end).anaWins = ana_win(lib);

%For ITR - uses only trajectory
analysisParams(end+1) = analysisParams(strcmpi({analysisParams.name}, 'tuning'));
analysisParams(end).name = 'ITR';
% analysisParams(end).useUnits = 'all';
[~, lib] = ismember({'trajectory'}, {ana_win.name});
analysisParams(end).anaWins = ana_win(lib);

%For context tuning, use visual and cue periods only
analysisParams(end+1) = analysisParams(strcmpi({analysisParams.name}, 'tuning'));
analysisParams(end).name = 'ContextTuning';
% analysisParams(end).useUnits = 'all';
[~, lib] = ismember({'baseline', 'visual', 'cue'}, {ana_win.name});
analysisParams(end).anaWins = ana_win(lib);
analysisParams(end).classifyTarg = 'combinedRegionRule';
%% Return only the expected analysisParams struct
if nargin > 0 && any(strcmpi({analysisParams.name}, varargin{1}))
    analysisParams = analysisParams(strcmpi({analysisParams.name}, varargin{1}));
else
    analysisParams = analysisParams(strcmpi({analysisParams.name}, 'default'));
end