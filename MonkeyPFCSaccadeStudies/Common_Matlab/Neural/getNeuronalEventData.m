function nevDat = getNeuronalEventData(dataFPrefix, filepath, varargin)
%function nevDat = getNeuronalEventData(dataFPrefix, filepath)
%function nevDat = getNeuronalEventData(dataFPrefix, filepath, sorter)
%where sorter is 'plx', 'mksort', 'guill', or 'nsn'.
%This function's input arguments are (1) the data file prefix, without its
%extension, and (2) the path to find the data files. It searches the path
%for .nev files matching the file name (fallback on .nevmat).
%If the optional argument (3) is passed then it will also search for
%spike-sorter outputs that match the expected format based on the type of
%spike-sorter indicated by the argument.
%If the spike-sorter output is found then it will be merged into the
%output.
%The output nevDat is a data structure as follows:
%TODO: specify the data structure.
%.timeUnits
%.words
%.spikeTimes.invalid
%.contNeuroInfo.invalid
%.spikeTimes(ee).unit
%.timeUnits
%.trialID

%% Find the datafile
completeFilePath = fullfile(filepath, [dataFPrefix '.nev']);


clearvars -global nevDat%for some reason nevDat is a global !?

if exist(completeFilePath, 'file')
    load('-mat', completeFilePath);
    %nevDat.
    %analogEvents contNeuroInfo converted2ms fileID fileStartTime filename
    %spikeTimes trialID words
else %TODO: Load the file from original format.
%   fp_getDataFromNEV
%       fp_getSpikeTimesFromNEV
%       fp_getAnalogEventsFromNEV
%       fp_getTrialID; useful if cerebus got trial id via dio during recording
%       fp_getFileID
%   fp_getDataFromNS
%       fp_getContinuousNeuroDataFromNEV
end

%% Note units are in sec
%fp_second2Millisecond
nevDat.timeUnits = 'sec';