function cDat = checkBaselineFR(cDat, varargin)
% nevDat = checkBaselineFR(cDat)
% nevDat = checkBaselineFR(cDat, params)
% params is a struct with fields
% .baselineFRLims [0 Inf]
% .rejectAnyBlock [false] - If true, reject a unit if it fails in any
% single block of trials.
% .critFlip - ['none'] use all samples, or if critFlip is any of the
% cDat.flipNames, then only consider data up to that point when calculating
% firing rate.
% Doesn't really check 'baseline', just calculates average FR.
% TODO: Only use pre-X rasters.
global DEF

params = varg2params(varargin, struct(...
    'baselineFRLims', [0 Inf],...
    'rejectAnyBlock', false,...
    'critFlip', 'none'),...
    {'baselineFRLims', 'rejectAnyBlock', 'critFlip'});


nUnits = length(cDat.invalidUnits);

% If we have critFlip, and we know flipNames, then only look inbetween the
% first flip screen and the critFlip flipScreen
if ~strcmpi(params.critFlip, 'none')...
        && isfield(cDat, 'flipNames') && ~isempty(cDat.flipNames)...
        && any(strcmpi(cDat.flipNames, params.critFlip))...
        && find(strcmpi(cDat.flipNames, params.critFlip), 1, 'first') > 1
    cf_ix = find(strcmpi(cDat.flipNames, params.critFlip), 1, 'first');
    unit_fired = zeros(1, nUnits);
    nSamps = 0;
    for tr_ix = 1:length(cDat.trial)
        tr = cDat.trial(tr_ix);
        if length(tr.flipScreen) > cf_ix
            t1 = ceil(tr.flipScreen(1));
            t2 = floor(tr.flipScreen(cf_ix));
            tr_unit_fired = full(sum(tr.raster(t1:t2, :), 1));
            unit_fired = unit_fired + tr_unit_fired;
            nSamps = nSamps + t2 - t1;
        end
    end
    firingRates = 1000.0*unit_fired./nSamps;
else % Use all samples
    raster = full(cat(1, cDat.trial(1:end-1).raster));
    [nSamples, ~] = size(raster);
    firingRates = sum(raster).*1000./nSamples;
end
frHighBool = firingRates' > params.baselineFRLims(2);
frLowBool = firingRates' < params.baselineFRLims(1);
hbBool = cDat.invalidUnits == DEF.unitState.valid & (frHighBool | frLowBool);
cDat.invalidUnits(hbBool) = DEF.unitState.highBaseline;
if any(firingRates') > 100
    fprintf('Found a unit with baseline > 100.\n');
end

% If cDat has a block-field, do the check per-block and eliminate any units
% that fail for any single block.
if params.rejectAnyBlock && isfield(cDat.trial, 'cueTargRuleGroup')
    block_id = [cDat.trial.cueTargRuleGroup];
    block_uq = unique(block_id);
    nBlocks = length(block_uq);
    block_firingRates = nan(nUnits, nBlocks);
    for block_ix = 1:length(block_uq)
        block_bool = block_id == block_uq(block_ix);
        raster = full(cat(1, cDat.trial(block_bool).raster));
        nSamples = size(raster, 1);
        block_firingRates(:, block_ix) = sum(raster).*1000./nSamples;
    end
    frHighBool = any(block_firingRates > params.baselineFRLims(2), 2);
    frLowBool = any(block_firingRates < params.baselineFRLims(1), 2);
    hbBool = cDat.invalidUnits == DEF.unitState.valid & (frHighBool | frLowBool);
    cDat.invalidUnits(hbBool) = DEF.unitState.highBaseline;
end