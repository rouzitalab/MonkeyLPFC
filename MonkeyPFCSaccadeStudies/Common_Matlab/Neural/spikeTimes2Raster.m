function nevDat = spikeTimes2Raster(nevDat)
%nevDat = spikeTimes2Raster(nevDat)
%Convert spikeTimes to a sparse nevDat.raster
%[times x electrodes*units]
%Also saves nevDat.rasterTimes (a timeStamp vector)

timeFac = 1;
if isfield(nevDat, 'timeUnits') && strcmpi(nevDat.timeUnits(1), 's')
    timeFac = 1000;
end
lastFileTimeMS  = int32(timeFac*nevDat.words(end,2));

invalidUnits = [nevDat.spikeTimes.invalid]';
[nElecs, nUnits] = size(invalidUnits);

if ~isempty(nevDat.contNeuroInfo)
    badElec = [nevDat.contNeuroInfo.invalid]~=0;
    invalidUnits(badElec,:) = -1;
end

goodUnits = invalidUnits==0;
raster = false(lastFileTimeMS, nElecs*nUnits);
for ee = 1:nElecs
    gdUnitIx = find(goodUnits(ee,:));
    for uu = 1:length(gdUnitIx)
        spikeTimes = nevDat.spikeTimes(ee).unit{gdUnitIx(uu)};
        spikeTimes = int32(timeFac*spikeTimes);
        spikeTimes(spikeTimes==0) = []; %Sometimes the first spike is at t=0.
        raster(spikeTimes, (gdUnitIx(uu)-1)*nElecs + ee) = true;
    end
end
raster = sparse(raster);
    
nevDat.raster = raster;
nevDat.invalidUnits = invalidUnits;
nevDat.lastFileTimeMS = lastFileTimeMS;
nevDat = rmfield(nevDat, {'spikeTimes'});