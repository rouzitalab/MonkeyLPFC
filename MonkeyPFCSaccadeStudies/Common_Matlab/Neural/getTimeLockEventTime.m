function ev_times = getTimeLockEventTime(cDat, varargin)
%cDat = getTrialEventTime(cDat, varargin)
params = varg2params(varargin,...
    struct(...
    'timeLockEvent', 'fixationOffset'));

% This will be used below to make sure each trial has enough data to
% support winEdges
nTrials = length(cDat.trial);
ev_times = nan(1, nTrials);
if strcmpi(params.timeLockEvent, 'saccadeOnset')
    ev_times =[cDat.trial.sacStartTime];
elseif sum(strcmpi(cDat.flipNames, params.timeLockEvent)) == 1
    flip_ix = find(strcmpi(cDat.flipNames, params.timeLockEvent));
    for tr_ix = 1:nTrials
        if length(cDat.trial(tr_ix).flipScreen) >= flip_ix
            ev_times(tr_ix) = cDat.trial(tr_ix).flipScreen(flip_ix); %#ok<*FNDSB>
        end
    end
end