function mergeNevPtbIntoMKSort(nevfull, ptbfull, mkdir)
%function mergeNevPtbIntoMKSort(nevfull, ptbfull, mkdir)
%nevfull is fullfile path to the .nev file.
%ptbfull is fullfile path to the ptb file.
%mkdir is the mksort directory.
% This function loads data from the nev file and from the ptbfile and
% merges some of its information into the mksort files. This is used only
% because mksort has the ability to display this information within it.
% Other spike sorters cannot store the same information in their output. It
% should be expected that the output of spike sorters will not have trial
% information (e.g., condition, synchronization, etc.)

nevFileInfo = getNevFileInfo(nevfull);
nTrials = length(nevFileInfo.trialInfo);
trialStarts = 1000*[nevFileInfo.trialInfo.startTime];
trialEnds = 1000*[nevFileInfo.trialInfo.stopTime];
if length(trialEnds) < nTrials || isnan(trialEnds(end))
    lti = nevFileInfo.trialInfo(end);
    last_event = max([lti.flipScreen'; lti.leverDown; lti.leverRelease]);
    trialEnds(nTrials) = 1000*last_event;
end

ptb = load(ptbfull, '-mat');

%TODO: conditions can be from .class or .targetChoice per-trial.
%       Does that ever happen?
if isfield(ptb.eData.trial, 'class')
    conditions = [ptb.eData.trial.class];
elseif isfield(ptb.eData.trial, 'targetChoice')
    conditions = [ptb.eData.trial.targetChoice];
end

wfns = dir(fullfile(mkdir, 'waveforms_*.mat'));
fprintf('Merging nevFileInfo into mksort waveforms located in %s\n', mkdir);
for wf = 1:length(wfns)
    %For each waveform, 
    %waveforms.trialInfo.condition = array of tasks
    %waveforms.trialInfo.trialStartTimes = array of start times
    %waveforms.trialInfo.trialEndTimes = array of end times
    %waveforms.trialInfo.trial(spike_indices) = trial_id;
    loadvar = load(fullfile(mkdir, wfns(wf).name));
    waveforms = loadvar.waveforms;
    waveforms.trialInfo.condition = conditions;
    waveforms.trialInfo.trialStartTimes = trialStarts;
    waveforms.trialInfo.trialEndTimes = trialStarts + trialEnds;
    waveforms.trialInfo.trial = nan(1, length(waveforms.spikeTimes));
    for tr = 1:nTrials
        waveforms.trialInfo.trial(waveforms.spikeTimes >= waveforms.trialInfo.trialStartTimes(tr) &...
            waveforms.spikeTimes <= waveforms.trialInfo.trialEndTimes(tr)) = tr;
    end
    save(fullfile(mkdir, wfns(wf).name), 'waveforms');
    fprintf('.');
end
fprintf('Done.\n');