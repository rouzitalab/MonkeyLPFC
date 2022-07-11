function raster = getEpochRaster(cDat, anaWin)

cDat = trimTrials(cDat, anaWin);
tVec = anaWin.winEdges(1):anaWin.winEdges(2);
raster = nan(length(cDat.trial), sum(~cDat.invalidUnits), length(tVec));
for tr_ix = 1:length(cDat.trial)
    rBool = cDat.trial(tr_ix).timeBool;  % this trial's in-window samples
    if any(rBool)
        trTVec = cDat.trial(tr_ix).tVec(rBool);  % the times of this trial's in-window samples
        lBool = tVec >= trTVec(1) & tVec <= trTVec(end);  % the samples in the common raster
        raster(tr_ix, :, lBool) = cDat.trial(tr_ix).raster(rBool, :)';
    end
end