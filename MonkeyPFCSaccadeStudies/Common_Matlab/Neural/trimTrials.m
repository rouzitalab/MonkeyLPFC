function cDat = trimTrials(cDat, varargin)
%cDat = trimTrials(cDat)
%cDat = trimTrials(cDat, params)
%params has fields
% .timeLockEvent - the name of the event 
% .winEdges [-2000 2000]
% .avoidEvent [nan]
% .avoidWin [Inf -Inf]
% .commonWin [false]

params = varg2params(varargin,...
    struct('timeLockEvent', 'saccade onset',...
    'winEdges', [-2000 2000],...
    'avoidEvent', nan,...
    'avoidWin', [Inf -Inf],...
    'commonWin', false),...
    {'timeLockEvent' 'winEdges' 'avoidEvent' 'avoidWin' 'commonWin'});

timeFac = 1;
if strcmpi(cDat.timeUnits(1), 's')
    timeFac = 1000;
end

%% Get timeLockEvent time
evTimes = getTimeLockEventTime(cDat, struct('timeLockEvent', params.timeLockEvent));
evTimes = round(evTimes);
evTimes = num2cell(evTimes);
[cDat.trial.eventTime] = evTimes{:};
clear evTimes

%% Get the time of the event to avoid
if isfield(cDat.trial, 'avoidTime')
    cDat.trial = rmfield(cDat.trial, 'avoidTime');
end
cDat.trial(1).avoidTime = [];
if ~isnan(params.avoidEvent)
    evTimes = getTimeLockEventTime(cDat, struct('timeLockEvent', params.avoidEvent));
    evTimes = round(evTimes);
    evTimes = num2cell(evTimes);
    [cDat.trial.avoidTime] = evTimes{:};
    clear evTimes
end

%% Get a tVec and identify good samples for each trial.

% %fieldnames containing times that will be re-zero'd.
% fns = fieldnames(cDat.trial);
% fns = fns(ismember(fns,...
%     {'stopTime', 'flipScreen', 'sacStartTime', 'eventTime', 'leverRelease', 'avoidTime'}));

for tr_ix = 1:length(cDat.trial)
    tr = cDat.trial(tr_ix);
    evTVec = 1:size(tr.raster, 1);
    evTVec = evTVec - round(tr.eventTime);
    timeBool = evTVec >= params.winEdges(1) & evTVec < params.winEdges(2);
    
    % set timeBool = false where it overlaps with avoidWindow
    if ischar(params.avoidEvent) && ~isempty(params.avoidEvent) && ~any(isnan(params.avoidWin))
        avTVec = 1:size(tr.raster, 1);
        avTVec = avTVec - tr.avoidTime;
        avTimeBool = avTVec >= params.avoidWin(1) & avTVec < params.avoidWin(2);
        timeBool = timeBool & ~avTimeBool;
    end
    
    cDat.trial(tr_ix).timeBool = timeBool;
    cDat.trial(tr_ix).tVec = evTVec;
end

%% Use a common minimum window, if desired.
if params.commonWin
    %Identify the limits for each trial around the timeLockEvent
    tLims = nan(length(cDat.trial), 2);
    for tr_ix = 1:length(cDat.trial)
        trTVec = cDat.trial(tr_ix).tVec(cDat.trial(tr_ix).timeBool);
        tLims(tr_ix, :) = [min(trTVec) max(trTVec)];
    end
    
    %Mark samples outside the most restrictive limits as bad.
    for tr_ix = 1:length(cDat.trial)
        cDat.trial(tr_ix).timeBool = cDat.trial(tr_ix).tVec >= max(tLims(:,1)) & cDat.trial(tr_ix).tVec <= min(tLims(:, 2));
    end
end