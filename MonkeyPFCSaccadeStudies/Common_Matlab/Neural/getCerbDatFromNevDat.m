function cDat = getCerbDatFromNevDat(nevDat, DIO)
%function cDat = getCerbDatFromNevDat(nevDat)

global lim %DEF

%Create the output variable.
cDat = struct(...
    'trial', [],...
    'invalidUnits', nevDat.invalidUnits,...
    'fileStopTime', [],...
    'fileAlignTime', [],...
    'cerebusSamplingRate', 30000,...
    'lastValidTrialTime', [],...
    'timeUnits', nevDat.timeUnits);

%Time (in s) the file stopped.
stopRecIx = find(nevDat.words(:,1) == DIO.stopRecordingWord);
if ~isempty(stopRecIx)
    cDat.fileStopTime = nevDat.words(stopRecIx,2);
end

[nElecs, nUnits] = size(cDat.invalidUnits);

%Time (in s) each trial started and stopped.
allTrialStartIx = find(nevDat.words(:,1)==DIO.startTrialWord);
allTrialStopIx = find(nevDat.words(:,1)==DIO.stopTrialWord);
allTrialStartIx(allTrialStartIx>allTrialStopIx(end)) = []; %starts after last stop
allTrialStopIx(allTrialStopIx<allTrialStartIx(1)) = []; %stops before first start
allStartTimes = nevDat.words(allTrialStartIx,2);
allStopTimes = nevDat.words(allTrialStopIx,2);
%end of last trial in s.
cDat.lastValidTrialTime = nevDat.words(allTrialStopIx(end),2);

%Event times.
allLeverDownTimes    = nevDat.words(nevDat.words(:,1)==DIO.leverDown,2);
allLeverReleaseTimes = nevDat.words(nevDat.words(:,1)==DIO.leverRelease,2);
allFlipScreenTimes   = nevDat.words(nevDat.words(:,1)==DIO.flipScreen,2);

timeFac = 1;
if strcmpi(cDat.timeUnits(1),'s')
    timeFac = 1000;
end

nTrials = length(allTrialStopIx);%Can check numTrials > lim.maxNumTrials

%Create the struct array for trials.
cDat.trial = struct(...
    'startTime', num2cell(allStartTimes),...
    'stopTime', num2cell(allStopTimes),...
    'ID', num2cell(nevDat.trialID(1:nTrials,:),2),...
    'invalid', num2cell(int16((allStopTimes-allStartTimes)>lim.maxTrialDuration)),...
    'alignTime', num2cell(allStartTimes),...
    'leverDown', [],... cannot automate because not 1 per trial
    'leverRelease', [],...
    'flipScreen', [],...
    'raster', []); %trials are of different length so raster will change.



fprintf('Create trial-structure for neurodata: n = %d trials; 50 trials per dot.\n',nTrials);
for tt=1:nTrials
    if ~mod(tt, 1000)
        fprintf('\n');
    end
    if ~mod(tt,50)
        fprintf('.');
    end
    
    thisLDT = allLeverDownTimes(allLeverDownTimes>cDat.trial(tt).startTime & allLeverDownTimes<cDat.trial(tt).stopTime);
    if ~isempty(thisLDT)
        cDat.trial(tt).leverDown = thisLDT(1) - cDat.trial(tt).alignTime;
    end
    
    thisLRT = allLeverReleaseTimes(allLeverReleaseTimes>cDat.trial(tt).startTime & allLeverReleaseTimes<cDat.trial(tt).stopTime);
    if ~isempty(thisLRT)
        cDat.trial(tt).leverRelease = thisLRT(1) - cDat.trial(tt).alignTime;
    end
    
    thisFlip = allFlipScreenTimes(allFlipScreenTimes>cDat.trial(tt).startTime & allFlipScreenTimes<cDat.trial(tt).stopTime);
    if ~isempty(thisFlip)
        cDat.trial(tt).flipScreen = thisFlip - cDat.trial(tt).alignTime;
    end
    
    trStart = floor(timeFac*cDat.trial(tt).startTime);
    trStop = ceil(timeFac*cDat.trial(tt).stopTime);
    cDat.trial(tt).raster = nevDat.raster(trStart:trStop,:);
end
fprintf('Done.\n');