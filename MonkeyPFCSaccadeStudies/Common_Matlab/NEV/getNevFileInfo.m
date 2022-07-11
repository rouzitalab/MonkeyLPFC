function fileInfo = getNevFileInfo(nevFile)
%function fileID = getFileID(nevFile)
%return the fileID of a nevFile.
global paths
addpath(genpath(fullfile(paths.npmk)));

NEV = openNEV(nevFile, 'noread', 'nosave');

fileInfo = struct('fileID', 'none', 'trialInfo', [],...
    'SampleFreq', NEV.MetaTags.SampleRes,...
    'DataDuration', NEV.MetaTags.DataDuration,...
    'DataDurationSec', NEV.MetaTags.DataDurationSec,...
    'DateTimeRaw', NEV.MetaTags.DateTimeRaw);

wordMatrix = zeros(length(NEV.Data.SerialDigitalIO.UnparsedData), 2, 'double');
if ~isempty(wordMatrix)
    wordMatrix(:,1) = NEV.Data.SerialDigitalIO.UnparsedData;
    wordMatrix(:,2) = NEV.Data.SerialDigitalIO.TimeStampSec';
end
%[DIO, sentenceID] = get_default_DIO_sentenceID();
%[sentences, events] = getInfoFromWords(wordMatrix, 'DIO', DIO, 'sentenceID', sentenceID);
[sentences, events] = getInfoFromWords(wordMatrix);

fid_sentence = sentences(strcmpi({sentences.type}, 'fileID'));
if ~isempty(fid_sentence)
    fileInfo.fileID = char(fid_sentence.value');
elseif isfield(NEV, 'MetaTags') && isfield(NEV.MetaTags, 'Filename')
    fileInfo.fileID = NEV.MetaTags.Filename;
end
    

%Maybe we have a stopRecordingWord. Its value is probably the same as
%NEV.MetaTags.DataDurationSec.
stop_recording_bool = strcmpi({events.type}, 'stopRecordingWord');
if any(stop_recording_bool)
    fileInfo.stopRecording = events(stop_recording_bool).time;
end
events(stop_recording_bool) = [];

sentences = sentences(strcmpi({sentences.type}, 'trialID'));
nTrials = length(sentences);
trialInfo = repmat(...
    struct(...
    'datenum', [],...
    'idx', [],...
    'flipScreen', [],...
    'startTime', [],...
    'stopTime', [],...
    'leverDown', [],...
    'leverRelease', [],...
    'startMyStim', [],...
    'fixating', [],...
    'breakFixation', []),...
    nTrials, 1);

event_start = find(strcmpi({events.type}, 'startTrialWord'));
event_stop = find(strcmpi({events.type}, 'stopTrialWord'));
for tt = 1:nTrials
    V = double(sentences(tt).value)';
    if V(1) < 1980
        V(1) = V(1) + 2000;
    end
    trialInfo(tt).datenum = datenum(V(1:6));
    trialInfo(tt).idx = V(7)*100 + V(8);
    
    th_si = event_start(tt); %this start index
    th_st = events(th_si).time; %this start time
    trialInfo(tt).startTime = th_st; %seconds since ???
    
    th_ei = find(event_stop > th_si, 1, 'first');
    if ~isempty(th_ei)
        th_ei = event_stop(th_ei);
        trialInfo(tt).stopTime = events(th_ei).time - th_st;
    else
        th_ei = length(events)+1;
        trialInfo(tt).stopTime = nan;
    end
    
    trial_events = events(th_si+1:th_ei-1);
    
    %Can sometimes end up with multiple startTrialWord events.
    uhoh_ix = find(strcmpi({trial_events.type}, 'startTrialWord'), 1, 'First');
    if ~isempty(uhoh_ix)
        trial_events = trial_events(1:uhoh_ix-1);
    end
    
    %Can sometimes end up with events with no type
    trial_events = trial_events(~cellfun(@isempty,{trial_events.type}));
    
    if length(trial_events) > 0
        fns = unique({trial_events.type});
        for fn = fns
            trialInfo(tt).(fn{1}) = [trial_events(strcmpi({trial_events.type}, fn)).time] - th_st;
        end
    end
end
    
fileInfo.trialInfo = trialInfo;

fprintf('Found %s with fileID %s and %i trials.\n', nevFile, fileInfo.fileID, length(fileInfo.trialInfo));