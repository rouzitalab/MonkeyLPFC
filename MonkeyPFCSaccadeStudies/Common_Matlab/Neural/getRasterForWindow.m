function raster = getRasterForWindow(cDat, ana_win)
%function raster = getRasterForWindow(cDat, ana_win)

% Identify eventTime, avoidTime, a tVec, and timeBool for each trial.
cDat = trimTrials(cDat, ana_win);

%Plug each trial's raster into a common raster.
tVec = ana_win.winEdges(1):ana_win.winEdges(2);
raster = nan(length(cDat.trial), length(tVec));
for tr_ix = 1:length(cDat.trial)
    rBool = cDat.trial(tr_ix).timeBool;  % this trial's in-window samples
    trTVec = cDat.trial(tr_ix).tVec(rBool);  % the times of this trial's in-window samples
    lBool = tVec >= trTVec(1) & tVec <= trTVec(end);  % the samples in the common raster
    raster(tr_ix, lBool) = cDat.trial(tr_ix).raster(rBool, unit_ix);
end