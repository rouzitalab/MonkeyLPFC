function ptb = getPTB(fullptbfn)
%function ptb = getPTB(filename)
% getPTB: opens a '.mat'-file created with the PsychophysicsToolbox
% (PTB).
% This attempts to normalize some of the things we need that differ across
% file types. See checkPTBCharaceteristics.m

%% Load ptbmat
ptb = load(fullptbfn, '-mat');

if ~isfield(ptb, 'DIO')  % Unnecessary they all have DIO
    [ptb.DIO, ptb.sentenceID] = get_default_DIO_sentenceID();
end

% if isfield(ptb, 'responseFM')
%     n_trials = length(ptb.responseFM);
%     trial_times = nan(n_trials,1);
%     for trial_ix = 1:n_trials
%         trial_times(trial_ix) = datenum(ptb.responseFM{trial_ix}.trialDate);
%     end
%     eData.trialDateNum = trial_times;
% end

%% Normalize some fields.

%eData.
%startEyeRecordings -> firstEyeSample, but there is no way to sync with the
%timer on the eyetracker! I removed these files.
if isfield(ptb, 'CALADJ')
    ptb.eData.eyePosCalibAdjust = ptb.CALADJ;
    ptb = rmfield(ptb, 'CALADJ');
end
%eData.trial
ptb.eData.trial = renamefield(ptb.eData.trial, 'trialStartTime', 'startTime');
ptb.eData.trial = renamefield(ptb.eData.trial, 'trialStopTime', 'stopTime');
if ~isfield(ptb.eData.trial(end), 'stopTime') || isempty(ptb.eData.trial(end).stopTime)
    ptb.eData.trial(end).stopTime = nan; %Fix lack of stopTime on last trial.
end
if (ptb.eData.trial(1).stopTime - ptb.eData.trial(1).startTime) >= 0
    newvals = num2cell([ptb.eData.trial.stopTime] - [ptb.eData.trial.startTime]);
    [ptb.eData.trial.stopTime] = newvals{:};  %Convert from absolute time to trial duration
end
ptb.eData.trial = renamefield(ptb.eData.trial, 'eyeSyncStartTime', 'eyeSyncTime'); %in eyetracker samples
if isfield(ptb.eData.trial, 'eyeSyncTime')
    newvals = num2cell(int32([ptb.eData.trial.eyeSyncTime]));
    [ptb.eData.trial.eyeSyncTime] = newvals{:};
end
ptb.eData.trial = renamefield(ptb.eData.trial, 'trialID', 'ID');
ptb.eData.trial = renamefield(ptb.eData.trial, 'timeofrelease', 'leverRelease');
ptb.eData.trial = renamefield(ptb.eData.trial, 'cueTime', 'cuePresentedTime');

if length(unique({ptb.eData.trial.expType})) > 1
    % In files with both M and A trial types,
    % M trials have targetChoice 1-8 (maybe 1-16) and class is empty.
    % A trials have targetChoice 0 or 1 and class is targetChoice+1.
    ptb.eData.trial = rmfield(ptb.eData.trial, 'class');
end
ptb.eData.trial = renamefield(ptb.eData.trial, 'targetChoice', 'class');

% trial.startTime is in seconds since some unknown time.
% all other values are in seconds since trial start.

%% words
% The ptb file has a wordsMatrix in each ptb.eData.trial(ix)
% Therein, there is the trialID sentence, with value = [y m d HH MM SS trial_hunds trial_mod_hund]
% There are event words. Typically startTrial, flipScreen (multiple),
% stopTrial.
% The startTrial's time is often the same as the trial's startTime (or trialStartTime) 
% The stopTrial's time is the same as the trial's startTime (or
% trialStartTime) + its stopTime (or trialStopTime).
% nTrials = length(ptb.eData.trial);
% my_trials = repmat(struct('startTime', [], 'stopTime', [], 'leverDownTime', [], ...
%     'leverReleaseTime', [], 'flipScreen', [], 'myStimTime', [],...
%     'fixating', [], 'breakFixationTime', [], 'trialID', []), nTrials, 1);
% for tr_ix = 1:nTrials
%     [sentences, events] = getInfoFromWords(ptb.eData.trial(tr_ix).words, 'DIO', ptb.DIO, 'sentenceID', ptb.sentenceID);
% end
end

%%
function my_struct = renamefield(my_struct, oldField, newField)
if isfield(my_struct, oldField) && ~isfield(my_struct, newField)
    [my_struct.(newField)] = my_struct.(oldField);
    my_struct = rmfield(my_struct, oldField);
end
end