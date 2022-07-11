function [cDat, fullInvalid] = triageUnits(cDat, varargin)
% [cDat, fullInvalid] = trimUnits(cDat)
% 
% [cDat, fullInvalid] = triageUnits(cDat, params)
% Where params is a structure with any of the fields
% .useUnits - ['all'], 'sorted', or 'merged'. 'all' keeps all threshold
% crossings and separately sorted units, 'sorted' keeps only sorted units,
% 'merged' puts all threshold crossings together, effectively undoing
% sorting.
% params is also passed to checkBaselineFR (see its help)

% Note that fullInvalid is already used to reduce .trial(tt).raster,
% .invalidUnits, and .unitChannelAssignment

global DEF
params = varg2params(varargin,...
    struct(...
    'baselineFRLims', [0 Inf],...
    'useUnits', 'all',...  % 'all' 'sorted' 'merged'
    'rejectAnyBlock', false,...
    'critFlip', 'none'),...
    {'baselineFRLims', 'maxAllowedBaselineFR', 'minAllowedBaselineFR', 'useUnits', 'rejectAnyBlock', 'critFlip'});

nUnits = length(cDat.unitChannelAssignment);
if ~isfield(cDat, 'invalidUnits')
    cDat.invalidUnits = ones(nUnits, 1)  * DEF.unitState.valid;
end

%% Keep all units, keep only sorted units, or merge all into thresh cross
if ~strcmpi(params.useUnits, 'all')
    uqChans = unique(cDat.unitChannelAssignment);
    nChans = length(uqChans);
    nUnits = length(cDat.unitChannelAssignment);

    if strcmpi(params.useUnits, 'merged')
        %Merge all units per channel to one thresh cross unit per channel
        unitMapBool = false(nUnits, nChans);
        for uu = 1:nUnits
            unitMapBool(uu, :) = ismember(uqChans, cDat.unitChannelAssignment(uu));
        end
        for tt = 1:length(cDat.trial)
            newRaster = false(size(cDat.trial(tt).raster, 1), nChans);
            for cc = 1:nChans
                newRaster(:, cc) = any(cDat.trial(tt).raster(:, unitMapBool(:, cc)), 2);
            end
            cDat.trial(tt).raster = sparse(newRaster);
        end
        % If merging units, the original cDat.unitChannelAssignment and
        % invalidateUnits no longer work.
        cDat.unitChannelAssignment = (1:nChans)';
        cDat.invalidUnits = ones(nChans, 1) * DEF.unitState.valid;
        
    elseif strcmpi(params.useUnits, 'sorted')
        %Don't use the first unit from each channel
        ubool = true(1, nUnits);
        for cc = 1:nChans
            ubool(find(cDat.unitChannelAssignment == uqChans(cc), 1, 'first')) = false;
        end
        cDat.invalidUnits(~ubool) = DEF.unitState.manualExclusion;
    end
end

%% Mark too slow or too fast firing unts as invalid.
cDat = checkBaselineFR(cDat, params);  % Uses DEF.unitState.highBaseline if above or below max/min

%% Actually remove marked units.
fullInvalid = cDat.invalidUnits;
unitBool = cDat.invalidUnits==0;
for tt = 1:length(cDat.trial)
    cDat.trial(tt).raster = cDat.trial(tt).raster(:, unitBool);
end
cDat.invalidUnits = cDat.invalidUnits(unitBool);
cDat.unitChannelAssignment = cDat.unitChannelAssignment(unitBool);
if isfield(cDat, 'sortQuality') && ~isempty(cDat.sortQuality)
    cDat.sortQuality = cDat.sortQuality(unitBool);
end
