function cDat = getCerbDat(this_sess, varargin)
%cDat = getCerbDat(this_sess)
%cDat = getCerbDat(this_sess, 'neur_pref', {'cDat', 'nevmat', 'plx', 'mksort', 'nev'})
%Loads neural data from data, trying to find instances of the data in the
%format specified in neur_pref (cell array in desc order).
%cDat is the final format (Chad's version). 
%nevmat is Florian's format.
%nev is the format in which it was recorded. Some online sorting already
%done.
%plx is the output of Plexon Offline Sorter.
%mksort is the output of Matt Kaufmann's spike sorter.
%The spike sorter outputs require the nev file to be present so that we may
%extract trial information.
%I have yet to implement nevmat, plx, or cDat


%% Defaults, Constants, and Parameters
global paths
params = struct(...
    'neur_pref', {{'cDat', 'nevmat', 'plx', 'mksort', 'nev'}});
params = varg2params(varargin, params, {'doSave', 'neur_pref'});
%% TODO: Figure out how big backPorch can be
backPorch = 300; % Keep up to XXX msec after the end of a trial stopTime.

% Prepare the output
cDat = struct(...
    'trial', [],...
    'fileStopTime', [],...
    'fileAlignTime', [],...
    'cerebusSamplingRate', [],...
    'lastValidTrialTime', [],...
    'unitChannelAssignment', [],...
    'sortQuality', [],...
    'timeUnits', 'msec');
%'invalidUnits', nevDat.invalidUnits,...

%% Determine which data format to import.
form_found = '';
for form_ix = 1:length(params.neur_pref)
    if strcmpi(params.neur_pref{form_ix}, 'cDat') && isempty(form_found)
        %TODO: Search for cDat
    elseif strcmpi(params.neur_pref{form_ix}, 'nevmat') && isempty(form_found)
        %TODO: Search for nevmat
    elseif strcmpi(params.neur_pref{form_ix}, 'nev') && isempty(form_found)
        nevfn = fullfile(paths.cdataRoot, this_sess.subject, this_sess.nevdir, this_sess.nevfname);
        if exist(nevfn, 'file') == 2
            form_found = 'nev';
        end
    elseif strcmpi(params.neur_pref{form_ix}, 'plx') && isempty(form_found)
        nevfn = fullfile(paths.cdataRoot, this_sess.subject, this_sess.nevdir, this_sess.nevfname);
        plxfn = fullfile(paths.preprocessed, 'plx', this_sess.subject, this_sess.nevdir, [this_sess.ptbfname(1:end-6) '.plx']);
    elseif strcmpi(params.neur_pref{form_ix}, 'mksort') && isempty(form_found)
        nevfn = fullfile(paths.cdataRoot, this_sess.subject, this_sess.nevdir, this_sess.nevfname);
        mksortdir = fullfile(paths.preprocessed, 'mk', this_sess.subject, ...
            this_sess.nevdir, this_sess.ptbfname(1:end-7));
        if exist(nevfn, 'file') == 2 && isdir(mksortdir)
            form_found = 'mksort';
        end
    end
end
if isempty(form_found)
    error('No data found matching preferred import list.')
end
clear form_ix

%% Load the trial information
%startTimes_msec, endTimes_msec, trial_ids
if strcmpi(form_found, 'cDat')
elseif strcmpi(form_found, 'nevmat')
    %TODO: nevmat->cDat
    %cDat = getCerbDatFromNevDat(nevDat, DIO)
else %nev, plx, mksort all get trial information from nev file.
    nevfn = fullfile(paths.cdataRoot, this_sess.subject, this_sess.nevdir, this_sess.nevfname);
    nevFileInfo = getNevFileInfo(nevfn);
    
    %per-file info
    cDat.cerebusSamplingRate = nevFileInfo.SampleFreq;
    
    %per-trial info
    startTimes_msec = [nevFileInfo.trialInfo.startTime]*1000;
    endTimes_msec = [nevFileInfo.trialInfo.stopTime]*1000;
    trial_ids = [nevFileInfo.trialInfo.idx];
    
    nTrials = length(nevFileInfo.trialInfo);
    if length(endTimes_msec) < nTrials
        endTimes_msec = [endTimes_msec nan];
    end
    %Create the struct array for trials.
    cDat.trial = struct(...
        'startTime', num2cell(startTimes_msec),...
        'stopTime', num2cell(endTimes_msec),...
        'ID', num2cell(int16(trial_ids)),...
        'invalid', num2cell(int16(false(1, nTrials))), ...%num2cell(int16((endTimes_msec-startTimes_msec)>lim.maxTrialDuration)),...
        'alignTime', num2cell(startTimes_msec),...
        'leverDown', [],... cannot automate because not 1 per trial
        'leverRelease', [],...
        'flipScreen', [],...
        'raster', []); %trials are of different length so raster will change.
    
    fns = {'flipScreen', 'leverDown', 'leverRelease', 'startMyStim', 'fixating', 'breakFixation'};
    for tr_ix = 1:nTrials
        for fn_ix = 1:length(fns)
            cDat.trial(tr_ix).(fns{fn_ix}) = 1000*nevFileInfo.trialInfo(tr_ix).(fns{fn_ix})';
        end
    end
    %event times are generally in msec since the start of the file.
    
    clear startTimes_msec endTimes_msec trial_ids tr_ix fns fn_ix
end

%% Load the unit information
%Output is cDat.unitChannelAssignment and raster
if strcmpi(form_found, 'cDat')
elseif strcmpi(form_found, 'nevmat')
    %TODO: nevmat->cDat
elseif strcmpi(form_found, 'nev')
    %nevfn = fullfile(paths.cdataRoot, this_sess.subject, this_sess.nevdir, this_sess.nevfname);
    [spikes, ~, units] = nev2MatSpikesOnly(nevfn);  % From mksort toolbox.
    %unit_idx = round(100*mod(units,1));
    cDat.unitChannelAssignment = floor(units);
    raster = false(ceil(1000*nevFileInfo.DataDurationSec), length(units));
    for unit_ix = 1:length(units)
        raster(round(1000*spikes(spikes(:,1) == units(unit_ix), 2)), unit_ix) = true;
    end
    
    clear spikes units unit_ix
elseif strcmpi(form_found, 'plx')
    %plxfn = fullfile(paths.preprocessed, 'plx', this_sess.subject, this_sess.nevdir, [this_sess.ptbfname(1:end-6) '.plx']);
    % TODO: Get unit information from plx
    %raster
    %cDat.unitChannelAssignment = ch_idx;
    %found_data = true;
elseif strcmpi(form_found, 'mksort')
    %mkdir = fullfile(paths.preprocessed, 'mk', this_sess.subject, this_sess.nevdir, this_sess.ptbfname(1:end-7));
    % Get unit information from mksort
    wfns = dir(fullfile(mksortdir, 'waveforms_*.mat'));
    wforms = cell(1,length(wfns));
    for wf = 1:length(wfns)
        loadVar = load(fullfile(mksortdir, wfns(wf).name));
        wforms{wf} = loadVar.waveforms;
    end
    nUnits = 0;
    for wf = 1:length(wforms)
        nUnits = nUnits + length(unique(wforms{wf}.units));
    end
    raster = false(ceil(wforms{wf}.trialInfo.trialEndTimes(end)), nUnits);
    ch_idx = nan(nUnits, 1);
    sort_qual = zeros(nUnits, 1);
    unit_ix = 0;
    for wf = 1:length(wforms)
        spike_unit_ids = unique(wforms{wf}.units);
        for uu = 1:length(spike_unit_ids)
            this_unit_bool = wforms{wf}.units == spike_unit_ids(uu);  % spikes with this unit label
            
            this_sptimes = round(wforms{wf}.spikeTimes(this_unit_bool));
            this_sptimes = this_sptimes(this_sptimes < size(raster, 1) & this_sptimes > 0);
            raster(this_sptimes, unit_ix+uu) = true;
            
            ch_idx(unit_ix+uu) = wf;
            if uu > 1
                sort_qual(unit_ix+uu) = wforms{wf}.ratings(1).ratings(uu-1);
            end
        end
        unit_ix = unit_ix + uu;
    end
    cDat.unitChannelAssignment = ch_idx;
    cDat.sortQuality = sort_qual;
    clear wfns wforms wf unit_ids uu this_sptimes ch_idx unit_ix sort_qual
end

%% Put trial information and unit information into cDat structure
if ~strcmpi(form_found, 'cDat')
    nTrials = length(cDat.trial);
    raster = sparse(raster);
    
    fprintf('Create trial-structure for cDat: n = %d trials.\n',nTrials);
    for tr_ix = 1:nTrials
        start_ix = round(cDat.trial(tr_ix).startTime);
        if ~isnan(cDat.trial(tr_ix).stopTime)
            %% TODO: Add some extra time, but don't go beyond the next trial startTime.
            stop_ix = round(cDat.trial(tr_ix).startTime + cDat.trial(tr_ix).stopTime + backPorch);
            if tr_ix < nTrials
                stop_ix = min(stop_ix, round(cDat.trial(tr_ix+1).startTime));
            else
                stop_ix = min(stop_ix, size(raster, 1));
            end
        elseif tr_ix < nTrials
            stop_ix = round(cDat.trial(tr_ix+1).startTime);
        else
            stop_ix = size(raster, 1);
        end
        cDat.trial(tr_ix).raster = raster(start_ix:stop_ix, :);
        if mod(tr_ix, 50)==0
            fprintf('.');
        end
    end
    fprintf('Done.\n');
end