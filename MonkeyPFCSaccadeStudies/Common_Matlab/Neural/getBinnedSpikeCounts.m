function cDat = getBinnedSpikeCounts(cDat, varargin)
%cDat = getBinnedSpikeCounts(cDat, varargin)
%Adds fields .binnedCount (nUnits x nTimes) and.tVec to each cDat.trial

params.binWidth = 20;
params.kernSD = 30;
params.doSqrt = true;
params.doSmooth = true;
params = varg2params(varargin, params, {'binWidth', 'kernSD', 'doSqrt', 'doSmooth'});

timeFac = 1;
if strcmpi(cDat.timeUnits(1), 's')
    timeFac = 1000;
end
nTrials = length(cDat.trial);

[cDat.trial.binnedCount] = deal(nan);
[cDat.trial.tVec] = deal(nan);

%% Counting spikes in each bin.
nGoodUnits = sum(cDat.invalidUnits==0);
fprintf('Calculating spikes for %i units in %i trials. (.=50 trials)',...
    nGoodUnits, nTrials);
tic;
for tt = 1:nTrials
    if mod(tt,1000)==1
        fprintf('\n');
    end
    if mod(tt,50)==0
        fprintf('.');
    end
    
    cTrial = cDat.trial(tt);
    
    %Get the trial-specific tvec that matches the .raster
    trTVec = (1:size(cTrial.raster, 1)) - timeFac*cTrial.eventTime;
    binEdges = 0:-params.binWidth:(min(trTVec) - params.binWidth + 1);
    binEdges = [fliplr(binEdges) params.binWidth:params.binWidth:max(trTVec)+params.binWidth-1];
    cTrial.tVec = binEdges(2:end) - params.binWidth/2; %Bin Centers
    nTimes = length(cTrial.tVec);
    
    %For each window, count the spikes
    cTrial.binnedCount = nan(nGoodUnits,nTimes);
    for xx = 1:nTimes
        xBool = trTVec >= binEdges(xx) & trTVec < binEdges(xx+1);
        if any(xBool) &&... %We have some samples
                trTVec(1) < binEdges(xx) &&... %We had some samples before this bin
                trTVec(end) > binEdges(xx+1) %We have some samples after this bin.
            cTrial.binnedCount(:, xx) = sum(cTrial.raster(xBool, :));
        end %otherwise bin is incomplete and spike count is unreliable.
    end
    
    if params.doSqrt
        cTrial.binnedCount = sqrt(cTrial.binnedCount);
    end
    if params.doSmooth
        cTrial.binnedCount = smoother(cTrial.binnedCount,...
            params.kernSD, params.binWidth, 'causal', true);
    end
    cDat.trial(tt) = cTrial;
end
telapsed = toc;
fprintf('Done in %i seconds.\n', round(telapsed));